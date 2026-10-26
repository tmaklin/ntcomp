// ntcomp: Sequencing data compression using SBWT and k-bounded matching statistics.
//
// Copyright 2025 Tommi Mäklin [tommi@maklin.fi].
//
// Copyrights in this project are retained by contributors. No copyright assignment
// is required to contribute to this project.
//
// Except as otherwise noted (below and/or in individual files), this
// project is licensed under the Apache License, Version 2.0
// <LICENSE-APACHE> or <http://www.apache.org/licenses/LICENSE-2.0> or
// the MIT license, <LICENSE-MIT> or <http://opensource.org/licenses/MIT>,
// at your option.
//
//! The archive's contig names. Assemblers and repositories name the contigs
//! of an assembly in a pattern (`NODE_12_length_5321_cov_8.21`,
//! `SRR5024312.3`, `contig00012 len=5321 ...`), so the names are coded as
//! tokens against the previous name, with the contig lengths and the input
//! file names, which the decoder knows before the names, as predictions.
//!
//! ## Stored names
//!
//! The section starts with its method: 0 for the names as LZMA (the
//! postcard list in [crate::codecs::lzma]), 1 for the token stream below.
//! The encoder writes both and keeps the shorter, so names without a
//! pattern cost no more than LZMA.
//!
//! The token stream codes the names file by file, all with the arithmetic
//! coder ([crate::codecs::arith]). Per file, an adaptive flag says the
//! order: the contigs' order or decreasing contig length (ties in the
//! contigs' order; SPAdes numbers its contigs in that order), whichever the
//! encoder finds cheaper. A name is split into tokens: runs of digits (at
//! most 18, as a number and its width, so leading zeros are kept), runs of
//! letters, and single other bytes. It is coded token by token against the
//! previous name of its file (a file's first name: the previous file's
//! first name). Per token position, one operation of eight (an adaptive
//! 3-bit tree, by the position and the operation at that position in the
//! previous name):
//! - MATCH: the previous name's token; END: no more tokens;
//! - DELTA: the previous name's number plus d ≥ 1 (adaptive Elias gamma), in
//!   the same width if that had leading zeros;
//! - LENGTH: the contig's length (with its runs of N and its edits);
//! - FILE: one of the first 16 tokens of the input file's name (a 4-bit
//!   tree);
//! - SELF: a number from up to 32 tokens back in the same name (a 5-bit
//!   tree), such as a length repeated in a long header;
//! - DIGITS: a new number: its width (a 5-bit tree) and its digits (4-bit
//!   trees by position);
//! - TEXT: new text: its length (adaptive Elias gamma) and bytes (8-bit
//!   trees by the previous byte).
//!
use std::ops::Range;

use super::{ChainCsr, Choices};
use super::snps::SnpEdits;
use crate::nruns::NRuns;
use crate::codecs::arith::{code_bit, code_tree, Coder, Decoder, Encoder, GammaModel, C0};
use crate::codecs::lzma;

type E = Box<dyn std::error::Error>;

/// The section's methods.
const METHOD_LZMA: u8 = 0;
const METHOD_TOKENS: u8 = 1;

/// Token positions with their own contexts; later ones share the last.
const N_POS: usize = 64;
/// The most tokens in a name and bytes in a token when decoding.
const MAX_TOKENS: usize = 1 << 16;
const MAX_TEXT: u64 = 1 << 16;
/// Tokens of the file name a FILE operation can name.
const N_FILE_TOKENS: usize = 16;
/// Tokens back a SELF operation can reach.
const SELF_BACK: usize = 32;

#[derive(Debug)]
pub struct NameDecodeError(pub String);

impl std::fmt::Display for NameDecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "contig names: {}", self.0)
    }
}

impl std::error::Error for NameDecodeError {}

#[derive(Clone, Debug, PartialEq)]
enum Token {
    Number { value: u64, width: u8 },
    Text(Vec<u8>),
}

impl Token {
    fn number(value: u64) -> Token {
        Token::Number { value, width: digits(value) }
    }
}

/// Digits of `value` without leading zeros.
fn digits(
    value: u64,
) -> u8 {
    value.checked_ilog10().map_or(1, |l| l as u8 + 1)
}

fn tokenize(
    name: &[u8],
) -> Vec<Token> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < name.len() {
        let mut j = i + 1;
        if name[i].is_ascii_digit() {
            while j < name.len() && name[j].is_ascii_digit() { j += 1; }
            if j - i <= 18 {
                let value = name[i..j].iter().fold(0_u64, |v, &d| v * 10 + (d - b'0') as u64);
                out.push(Token::Number { value, width: (j - i) as u8 });
            } else {
                out.push(Token::Text(name[i..j].to_vec()));
            }
        } else if name[i].is_ascii_alphabetic() {
            while j < name.len() && name[j].is_ascii_alphabetic() { j += 1; }
            out.push(Token::Text(name[i..j].to_vec()));
        } else {
            out.push(Token::Text(vec![name[i]]));
        }
        i = j;
    }
    out
}

fn render(
    tokens: &[Token],
) -> Vec<u8> {
    let mut out = Vec::new();
    for token in tokens {
        match token {
            Token::Number { value, width } => out.extend(format!("{:0w$}", value, w = *width as usize).bytes()),
            Token::Text(bytes) => out.extend_from_slice(bytes),
        }
    }
    out
}

const MATCH: u32 = 0;
const END: u32 = 1;
const DELTA: u32 = 2;
const LENGTH: u32 = 3;
const FILE: u32 = 4;
const SELF: u32 = 5;
const DIGITS: u32 = 6;
const TEXT: u32 = 7;
/// The operation context for a position the previous name does not have.
const NO_OP: usize = 8;

/// The adaptive models of the token stream.
#[derive(Clone)]
struct Models {
    order: u32,
    op: Vec<u32>,
    delta: Vec<GammaModel>,
    file: Vec<u32>,
    back: Vec<u32>,
    width: Vec<u32>,
    digit: Vec<u32>,
    text_len: Vec<GammaModel>,
    byte: Vec<u32>,
}

impl Models {
    fn new() -> Self {
        Models {
            order: C0,
            op: vec![C0; N_POS * (NO_OP + 1) * 8],
            delta: (0..N_POS).map(|_| GammaModel::new()).collect(),
            file: vec![C0; N_POS * 16],
            back: vec![C0; N_POS * 32],
            width: vec![C0; N_POS * 32],
            digit: vec![C0; N_POS * 20 * 16],
            text_len: (0..N_POS).map(|_| GammaModel::new()).collect(),
            byte: vec![C0; 256 * 256],
        }
    }
}

/// A name as the reference of the next: its tokens and their operations.
#[derive(Clone, Default)]
struct Coded {
    tokens: Vec<Token>,
    ops: Vec<u32>,
}

/// The operation coding `token` at `pos` of `name`, and its argument.
fn choose(
    token: &Token,
    name: &[Token],
    pos: usize,
    reference: Option<&Token>,
    length: u64,
    file: &[Token],
) -> (u32, u64) {
    if Some(token) == reference {
        return (MATCH, 0)
    }
    let Token::Number { value, width } = token else {
        return match file.iter().take(N_FILE_TOKENS).position(|f| f == token) {
            Some(j) => (FILE, j as u64),
            None => (TEXT, 0),
        }
    };
    if *value == length && *width == digits(length) {
        return (LENGTH, 0)
    }
    if let Some(back) = (1..=pos.min(SELF_BACK)).find(|&d| &name[pos - d] == token) {
        return (SELF, back as u64)
    }
    if let Some(Token::Number { value: r, width: rw }) = reference
        && value > r && value - r < 1 << 32
        && (width == rw || (*rw == digits(*r) && *width == digits(*value))) {
        return (DELTA, value - r)
    }
    if let Some(j) = file.iter().take(N_FILE_TOKENS).position(|f| f == token) {
        return (FILE, j as u64)
    }
    (DIGITS, 0)
}

/// Codes (or decodes, `name` None) a name against `reference`, the contig
/// `length` and the file name's tokens. None if the stream is corrupt.
fn code_name<C: Coder>(
    c: &mut C,
    m: &mut Models,
    name: Option<&[Token]>,
    reference: &Coded,
    length: u64,
    file: &[Token],
) -> Option<Coded> {
    let mut coded = Coded::default();
    for pos in 0..=MAX_TOKENS {
        if c.overrun() || pos == MAX_TOKENS {
            return None
        }
        let p = pos.min(N_POS - 1);
        let r = reference.tokens.get(pos);
        let (op, arg) = match name {
            None => (0, 0),
            Some(name) => match name.get(pos) {
                None => (END, 0),
                Some(token) => choose(token, name, pos, r, length, file),
            },
        };
        let previous = reference.ops.get(pos).map_or(NO_OP, |&o| o as usize);
        let op = code_tree(c, &mut m.op[(p * (NO_OP + 1) + previous) * 8..][..8], 3, op);
        let token = match op {
            END => return Some(coded),
            MATCH => r?.clone(),
            LENGTH => Token::number(length),
            DELTA => {
                let d = m.delta[p].code(c, arg);
                let Some(Token::Number { value, width }) = r else { return None };
                let v = value.checked_add(d)?;
                let width = if *width == digits(*value) { digits(v) } else { *width };
                if digits(v) > width || width > 18 { return None }
                Token::Number { value: v, width }
            },
            FILE => file.get(code_tree(c, &mut m.file[p * 16..][..16], 4, arg as u32) as usize)?.clone(),
            SELF => {
                let back = code_tree(c, &mut m.back[p * 32..][..32], 5, arg.saturating_sub(1) as u32) as usize + 1;
                coded.tokens.get(pos.checked_sub(back)?)?.clone()
            },
            DIGITS => {
                let truth = name.map(|n| render(&n[pos..pos + 1]));
                let width = code_tree(c, &mut m.width[p * 32..][..32], 5, truth.as_ref().map_or(0, |t| t.len() as u32)) as usize;
                if width == 0 || width > 18 { return None }
                let mut value = 0_u64;
                for i in 0..width {
                    let d = code_tree(c, &mut m.digit[(p * 20 + i) * 16..][..16], 4, truth.as_ref().map_or(0, |t| (t[i] - b'0') as u32));
                    if d > 9 { return None }
                    value = value * 10 + d as u64;
                }
                Token::Number { value, width: width as u8 }
            },
            _ => {
                let truth = name.map(|n| render(&n[pos..pos + 1]));
                let len = m.text_len[p].code(c, truth.as_ref().map_or(1, |t| t.len() as u64));
                if len > MAX_TEXT { return None }
                let mut bytes = Vec::with_capacity(len as usize);
                let mut previous = 0_usize;
                for i in 0..len as usize {
                    let b = code_tree(c, &mut m.byte[previous * 256..][..256], 8, truth.as_ref().map_or(0, |t| t[i] as u32));
                    bytes.push(b as u8);
                    previous = b as usize;
                }
                Token::Text(bytes)
            },
        };
        coded.tokens.push(token);
        coded.ops.push(op);
    }
    None
}

/// A coder that only adds up the cost of the bits, for choosing an order.
struct Cost(f64);

impl Coder for Cost {
    fn code(&mut self, bit: u32, p1: u32) -> u32 {
        let p = p1 as f64 / 65536.0;
        self.0 -= (if bit == 1 { p } else { 1.0 - p }).log2();
        bit
    }
}

/// The colours of `range` in decreasing length (ties in colour order) if
/// `by_length`, otherwise in colour order.
fn file_order(
    range: &Range<u32>,
    lengths: &[u64],
    by_length: bool,
) -> Vec<u32> {
    let mut order: Vec<u32> = range.clone().collect();
    if by_length {
        order.sort_by_key(|&c| std::cmp::Reverse(lengths[c as usize]));
    }
    order
}

/// The token stream of `names`, see the module docs.
fn encode_tokens(
    names: &[Vec<u8>],
    files: &[(Vec<u8>, Range<u32>)],
    lengths: &[u64],
) -> Vec<u8> {
    let mut encoder = Encoder::new();
    let mut m = Models::new();
    let mut first = Coded::default();
    let tokens: Vec<Vec<Token>> = names.iter().map(|n| tokenize(n)).collect();
    for (file, range) in files {
        let file = tokenize(file);
        // Code the file both ways with a copy of the models, keep the cheaper.
        let cost = |by_length: bool| {
            let (mut cost, mut m, mut reference) = (Cost(0.0), m.clone(), first.clone());
            for c in file_order(range, lengths, by_length) {
                reference = code_name(&mut cost, &mut m, Some(&tokens[c as usize]), &reference, lengths[c as usize], &file).expect("names encode");
            }
            cost.0
        };
        let by_length = range.len() > 1 && cost(true) < cost(false);
        code_bit(&mut encoder, &mut m.order, by_length as u32);
        let mut reference = first.clone();
        for (k, c) in file_order(range, lengths, by_length).into_iter().enumerate() {
            reference = code_name(&mut encoder, &mut m, Some(&tokens[c as usize]), &reference, lengths[c as usize], &file).expect("names encode");
            if k == 0 { first = reference.clone(); }
        }
    }
    encoder.finish()
}

/// Decodes the token stream of the names of `files`.
fn decode_tokens(
    bytes: &[u8],
    files: &[(Vec<u8>, Range<u32>)],
    lengths: &[u64],
) -> Result<Vec<Vec<u8>>, NameDecodeError> {
    let corrupt = || NameDecodeError("token stream truncated or corrupt".to_string());
    let mut decoder = Decoder::new(bytes);
    let mut m = Models::new();
    let mut first = Coded::default();
    let mut names: Vec<Vec<u8>> = vec![Vec::new(); lengths.len()];
    for (file, range) in files {
        let file = tokenize(file);
        let by_length = code_bit(&mut decoder, &mut m.order, 0) == 1;
        let mut reference = first.clone();
        for (k, c) in file_order(range, lengths, by_length).into_iter().enumerate() {
            reference = code_name(&mut decoder, &mut m, None, &reference, lengths[c as usize], &file).ok_or_else(corrupt)?;
            names[c as usize] = render(&reference.tokens);
            if k == 0 { first = reference.clone(); }
        }
    }
    if decoder.overrun() {
        return Err(corrupt())
    }
    Ok(names)
}

/// Every contig's length before it was stored (the names' LENGTH
/// predictions): the stored length (out of the source, see
/// [ChainCsr::sequence_len]) with its edits reverted and its runs of `N` put
/// back.
pub fn contig_lengths(
    graph: &ChainCsr,
    choices: &Choices,
    snp_edits: Option<&SnpEdits>,
    n_runs: Option<&NRuns>,
) -> Result<Vec<u64>, E> {
    (0..choices.streams.len() as u32).map(|color| {
        let stored = graph.sequence_len(choices, color)? as i64;
        let edits = snp_edits.map_or(Some(0), |s| s.length_change(color)).ok_or_else(|| NameDecodeError(format!("no edits of colour {}", color)))?;
        let runs: u64 = n_runs.map_or(Some(0), |r| r.runs.get(color as usize).map(|r| r.iter().map(|r| r.len as u64).sum())).ok_or_else(|| NameDecodeError(format!("no runs of colour {}", color)))?;
        let length = u64::try_from(stored + edits).map_err(|_| NameDecodeError(format!("negative length of colour {}", color)))?;
        Ok(length + runs)
    }).collect()
}

/// The files to code the names by: `files` if given, otherwise one unnamed
/// file of all `n` contigs. An error if they do not cover the contigs.
fn files_or_all(
    files: Option<&[(Vec<u8>, Range<u32>)]>,
    n: usize,
) -> Result<Vec<(Vec<u8>, Range<u32>)>, String> {
    let n = u32::try_from(n).map_err(|_| "more than 2^32 contigs")?;
    let files = files.map_or_else(|| vec![(Vec::new(), 0..n)], |f| f.to_vec());
    if files.last().map_or(0, |(_, r)| r.end) != n || files.windows(2).any(|w| w[0].1.end != w[1].1.start) || files.first().is_some_and(|(_, r)| r.start != 0) {
        return Err("the file ranges do not cover the contigs".to_string())
    }
    Ok(files)
}

/// The names section for `names` (one per contig), coded by the input
/// `files` (None: one) and the contigs' `lengths`, see the module docs.
pub fn encode(
    names: &[Vec<u8>],
    files: Option<&[(Vec<u8>, Range<u32>)]>,
    lengths: &[u64],
) -> Result<Vec<u8>, E> {
    if lengths.len() != names.len() {
        return Err(format!("{} contig lengths for {} names", lengths.len(), names.len()).into())
    }
    let files = files_or_all(files, names.len())?;
    // Both methods at once (LZMA takes longer).
    let (tokens, lzma) = rayon::join(
        || [vec![METHOD_TOKENS], encode_tokens(names, &files, lengths)].concat(),
        || -> Result<Vec<u8>, String> {
            let data = lzma::compress(&postcard::to_allocvec(names).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
            Ok([vec![METHOD_LZMA], postcard::to_allocvec(&data).map_err(|e| e.to_string())?].concat())
        });
    let lzma = lzma?;
    Ok(if tokens.len() <= lzma.len() { tokens } else { lzma })
}

/// Reads the names section written by [encode] with the same `files` and
/// `lengths`.
pub fn decode(
    bytes: &[u8],
    files: Option<&[(Vec<u8>, Range<u32>)]>,
    lengths: &[u64],
) -> Result<Vec<Vec<u8>>, E> {
    let files = files_or_all(files, lengths.len()).map_err(NameDecodeError)?;
    let names = match bytes.split_first() {
        Some((&METHOD_TOKENS, stream)) => decode_tokens(stream, &files, lengths)?,
        Some((&METHOD_LZMA, data)) => postcard::from_bytes(&lzma::decompress(&postcard::from_bytes(data)?)?)?,
        _ => return Err(Box::new(NameDecodeError("unknown method".to_string()))),
    };
    if names.len() != lengths.len() {
        return Err(Box::new(NameDecodeError(format!("{} names for {} contigs", names.len(), lengths.len()))))
    }
    Ok(names)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(names: &[&str], files: Option<&[(Vec<u8>, Range<u32>)]>, lengths: &[u64]) -> Vec<u8> {
        let names: Vec<Vec<u8>> = names.iter().map(|n| n.as_bytes().to_vec()).collect();
        let bytes = encode(&names, files, lengths).unwrap();
        assert_eq!(decode(&bytes, files, lengths).unwrap(), names);
        bytes
    }

    #[test]
    fn tokens_render_back() {
        for name in ["NODE_1_length_5321_cov_8.21", "contig00012 len=5321", "a", "", "x1234567890123456789012y", "ü 007 .-.", "0"] {
            let tokens = tokenize(name.as_bytes());
            assert_eq!(render(&tokens), name.as_bytes());
        }
        assert_eq!(tokenize(b"ab12.007"), vec![Token::Text(b"ab".to_vec()), Token::number(12), Token::Text(b".".to_vec()), Token::Number { value: 7, width: 3 }]);
    }

    #[test]
    fn patterned_names_use_the_tokens() {
        // SPAdes names, listed out of length order, and counter names.
        let lengths: Vec<u64> = vec![5000, 90000, 120, 70000, 300, 1000, 2000, 3000];
        let names = ["NODE_4_length_5000_cov_8.21_ID_7", "NODE_2_length_90000_cov_9.5_ID_3", "NODE_6_length_120_cov_1.0_ID_11",
                     "NODE_3_length_70000_cov_7.75_ID_5", "NODE_5_length_300_cov_2.25_ID_9",
                     "SRR5.1", "SRR5.2", "SRR5.003"];
        let files = vec![(b"spades".to_vec(), 0..5), (b"SRR5".to_vec(), 5..8)];
        let bytes = roundtrip(&names, Some(&files), &lengths);
        assert_eq!(bytes[0], METHOD_TOKENS);
        // Without files, and with no names.
        roundtrip(&names, None, &lengths);
        roundtrip(&[], None, &[]);
        roundtrip(&["", "", ""], Some(&[(b"f".to_vec(), 0..3)]), &[1, 2, 3]);
    }

    #[test]
    fn unpatterned_names_fall_back_to_lzma() {
        // Random identifiers: the tokens cost more than LZMA.
        let mut x = 88172645463325252_u64;
        let names: Vec<String> = (0..2000).map(|_| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; format!("{:016x}{:x}", x, x.rotate_left(17)) }).collect();
        let names: Vec<&str> = names.iter().map(|s| s.as_str()).collect();
        let lengths = vec![0_u64; names.len()];
        let bytes = roundtrip(&names, None, &lengths);
        assert_eq!(bytes[0], METHOD_LZMA);
    }

    #[test]
    fn corrupt_sections_fail() {
        let names = vec![b"NODE_1_length_5_cov_1.0".to_vec(); 3];
        let lengths = vec![5_u64; 3];
        let bytes = encode(&names, None, &lengths).unwrap();
        assert!(decode(&bytes[..bytes.len() / 2], None, &lengths).is_err());
        assert!(decode(&[9], None, &lengths).is_err());
        assert!(decode(&bytes, Some(&[(Vec::new(), 0..2)]), &lengths).is_err());
    }
}
