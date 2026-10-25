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
//! Runs of `N`, and of the other characters than A, C, G and T that a
//! normalized sequence can have ([OTHER]: IUPAC codes and gaps; nearly
//! always `N`). The index has none of them, so the archive cuts the runs
//! out of every contig before editing and parsing it, handles the pieces
//! between them as separate sequences ([segments]), and stores the runs to
//! put back after decoding ([NRuns]).
//!
//! ## Stored runs
//!
//! All colours are coded as one stream with the arithmetic coder
//! ([crate::codecs::arith]): per colour, the number of runs (adaptive Elias
//! gamma, by whether the previous colour had runs), then per run its
//! position in the contig without its runs, as the distance from the
//! previous run's position (adaptive Elias gamma), its length: a flag for
//! the previous run's length (assemblers write gaps of a fixed length),
//! otherwise adaptive Elias gamma; and a flag for a run of only `N`,
//! otherwise its characters, each as its index in [OTHER] (an adaptive
//! 4-bit tree).
//!
use std::ops::Range;

use crate::codecs::arith::{code_bit, code_tree, Coder, Decoder, Encoder, GammaModel, C0};
use crate::codecs::prealloc;

/// The characters other than A, C, G and T in a normalized sequence
/// ([crate::fastx::read_contigs]).
pub const OTHER: &[u8; 12] = b"N-RYSWKMBDHV";

/// A run of `len` characters other than A, C, G and T, before the base at
/// `pos` of the contig without its runs: `N`s if `bases` is empty,
/// otherwise `bases`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NRun {
    pub pos: u32,
    pub len: u32,
    pub bases: Vec<u8>,
}

/// Whether `base` is cut out (not A, C, G or T).
fn is_other(
    base: u8,
) -> bool {
    !matches!(base, b'A' | b'C' | b'G' | b'T')
}

#[derive(Debug)]
pub struct NRunDecodeError(pub String);

impl std::fmt::Display for NRunDecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "N runs: {}", self.0)
    }
}

impl std::error::Error for NRunDecodeError {}

/// `sequence` without its runs of characters other than A, C, G and T, and
/// the runs. None if a run or the sequence does not fit in u32, or a
/// character is not in [OTHER].
pub fn cut(
    sequence: &[u8],
) -> Option<(Vec<u8>, Vec<NRun>)> {
    if !sequence.iter().any(|&b| is_other(b)) {
        return Some((sequence.to_vec(), Vec::new()))
    }
    let mut out = Vec::with_capacity(sequence.len());
    let mut runs: Vec<NRun> = Vec::new();
    let mut i = 0;
    while i < sequence.len() {
        if !is_other(sequence[i]) {
            out.push(sequence[i]);
            i += 1;
            continue
        }
        let end = sequence[i..].iter().position(|&b| !is_other(b)).map_or(sequence.len(), |n| i + n);
        let run = &sequence[i..end];
        if !run.iter().all(|b| OTHER.contains(b)) {
            return None
        }
        let bases = if run.iter().all(|&b| b == b'N') { Vec::new() } else { run.to_vec() };
        runs.push(NRun { pos: out.len().try_into().ok()?, len: run.len().try_into().ok()?, bases });
        i = end;
    }
    Some((out, runs))
}

/// Puts `runs` back into `cut`, the sequence without them. None if they do
/// not fit it (not in order or past its end).
pub fn restore(
    cut: &[u8],
    runs: &[NRun],
) -> Option<Vec<u8>> {
    if runs.is_empty() {
        return Some(cut.to_vec())
    }
    let total = runs.iter().try_fold(cut.len(), |a, r| a.checked_add(r.len as usize))?;
    let mut out = Vec::with_capacity(total);
    let mut from = 0;
    for run in runs {
        let pos = run.pos as usize;
        if pos < from || pos > cut.len() {
            return None
        }
        out.extend_from_slice(&cut[from..pos]);
        if run.bases.is_empty() {
            out.resize(out.len() + run.len as usize, b'N');
        } else if run.bases.len() == run.len as usize {
            out.extend_from_slice(&run.bases);
        } else {
            return None
        }
        from = pos;
    }
    out.extend_from_slice(&cut[from..]);
    Some(out)
}

/// The pieces of a sequence of `len` bases without its `runs` that lie
/// between the runs: the nonempty ranges between the run positions.
pub fn segments(
    len: usize,
    runs: &[NRun],
) -> Vec<Range<usize>> {
    let mut out = Vec::with_capacity(runs.len() + 1);
    let mut start = 0;
    for run in runs {
        let pos = run.pos as usize;
        if pos > start { out.push(start..pos); }
        start = pos;
    }
    if len > start { out.push(start..len); }
    out
}

/// The runs of every colour, as stored in the archive.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NRuns {
    pub runs: Vec<Vec<NRun>>,
}

/// The adaptive models of the stream, see the module docs.
struct Models {
    count: [GammaModel; 2],
    gap: GammaModel,
    same_len: u32,
    len: GammaModel,
    only_n: u32,
    other: Vec<u32>,
}

/// Codes the runs of every colour (or decodes those of `n_colors` colours,
/// `runs` None). Decoding stops at an overrun or an invalid run.
fn code_all<C: Coder>(
    c: &mut C,
    n_colors: usize,
    runs: Option<&[Vec<NRun>]>,
) -> Result<Vec<Vec<NRun>>, NRunDecodeError> {
    let mut m = Models { count: [GammaModel::new(), GammaModel::new()], gap: GammaModel::new(), same_len: C0, len: GammaModel::new(), only_n: C0, other: vec![C0; 16] };
    let mut out: Vec<Vec<NRun>> = Vec::with_capacity(prealloc(n_colors));
    let mut had_runs = 0;
    let mut previous_len = 0_u64;
    for color in 0..n_colors {
        let truth = runs.map(|r| &r[color]);
        let n = m.count[had_runs].code(c, truth.map_or(1, |t| t.len() as u64 + 1)) - 1;
        if c.overrun() {
            return Err(NRunDecodeError(format!("truncated at colour {}", color)))
        }
        had_runs = (n > 0) as usize;
        let mut colour: Vec<NRun> = Vec::with_capacity(prealloc(n as usize));
        let mut pos = 0_u64;
        for i in 0..n as usize {
            let run = truth.map(|t| &t[i]);
            pos += m.gap.code(c, run.map_or(1, |r| r.pos as u64 - pos + 1)) - 1;
            let same = code_bit(c, &mut m.same_len, run.map_or(0, |r| (r.len as u64 == previous_len) as u32));
            let len = if same == 1 { previous_len } else { m.len.code(c, run.map_or(1, |r| r.len as u64)) };
            if c.overrun() || pos > u32::MAX as u64 || len == 0 || len > u32::MAX as u64 {
                return Err(NRunDecodeError(format!("invalid or truncated run of colour {}", color)))
            }
            previous_len = len;
            let only_n = code_bit(c, &mut m.only_n, run.map_or(0, |r| r.bases.is_empty() as u32));
            let mut bases = Vec::new();
            if only_n == 0 {
                for j in 0..len as usize {
                    let truth = run.map_or(0, |r| OTHER.iter().position(|&o| o == r.bases[j]).expect("a character of OTHER") as u32);
                    let symbol = code_tree(c, &mut m.other, 4, truth) as usize;
                    if c.overrun() || symbol >= OTHER.len() {
                        return Err(NRunDecodeError(format!("invalid or truncated run of colour {}", color)))
                    }
                    bases.push(OTHER[symbol]);
                }
            }
            colour.push(NRun { pos: pos as u32, len: len as u32, bases });
        }
        out.push(colour);
    }
    Ok(out)
}

impl NRuns {
    /// Number of colours.
    pub fn n_colors(&self) -> usize {
        self.runs.len()
    }

    /// Number of runs of all colours.
    pub fn n_runs(&self) -> usize {
        self.runs.iter().map(|r| r.len()).sum()
    }

    /// The runs coded with the arithmetic coder, see the module docs; the
    /// number of colours is stored elsewhere.
    pub fn to_bytes(
        &self,
    ) -> Vec<u8> {
        let mut encoder = Encoder::new();
        code_all(&mut encoder, self.runs.len(), Some(&self.runs)).expect("valid runs");
        encoder.finish()
    }

    /// Reads the runs of `n_colors` colours written by [NRuns::to_bytes].
    pub fn from_bytes(
        bytes: &[u8],
        n_colors: usize,
    ) -> Result<Self, NRunDecodeError> {
        let mut decoder = Decoder::new(bytes);
        let runs = code_all(&mut decoder, n_colors, None)?;
        if decoder.overrun() {
            return Err(NRunDecodeError("stream truncated".to_string()))
        }
        Ok(NRuns { runs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(pos: u32, len: u32) -> NRun {
        NRun { pos, len, bases: Vec::new() }
    }

    #[test]
    fn cut_and_restore() {
        for seq in [&b"ACGT"[..], b"", b"N", b"NNACGNNNTTN", b"ANA", b"NNNN", b"ACNNNNNNNNNNG", b"ACRYGNNSN-TK"] {
            let (cut, runs) = cut(seq).unwrap();
            assert!(cut.iter().all(|b| b"ACGT".contains(b)));
            assert_eq!(restore(&cut, &runs).unwrap(), seq);
            let pieces: Vec<u8> = segments(cut.len(), &runs).into_iter().flat_map(|r| cut[r].to_vec()).collect();
            assert_eq!(pieces, cut);
        }
        let (cut, runs) = cut(b"NNACGNNNTTN").unwrap();
        assert_eq!(cut, b"ACGTT");
        assert_eq!(runs, vec![run(0, 2), run(3, 3), run(5, 1)]);
        assert_eq!(segments(5, &runs), vec![0..3, 3..5]);
        // A run with other characters keeps them.
        let (_, runs) = super::cut(b"ACRNYGT").unwrap();
        assert_eq!(runs, vec![NRun { pos: 2, len: 3, bases: b"RNY".to_vec() }]);
        assert!(super::cut(b"ACXGT").is_none());
        assert!(restore(b"ACG", &[run(4, 1)]).is_none());
        assert!(restore(b"ACG", &[run(2, 1), run(1, 1)]).is_none());
        assert!(restore(b"ACG", &[NRun { pos: 1, len: 2, bases: b"R".to_vec() }]).is_none());
    }

    #[test]
    fn stream_roundtrip() {
        let runs = NRuns { runs: vec![
            vec![],
            vec![run(0, 250), run(1000, 250), NRun { pos: 1000, len: 3, bases: b"R-N".to_vec() }, run(1200, 1)],
            vec![run(5_000_000, 897), NRun { pos: 5_000_001, len: 1, bases: b"V".to_vec() }],
            vec![],
        ] };
        let bytes = runs.to_bytes();
        assert_eq!(NRuns::from_bytes(&bytes, 4).unwrap(), runs);
        assert!(NRuns::from_bytes(&bytes[..bytes.len() / 2], 4).is_err());
        // A corrupt colour count stops at the end of the stream.
        assert!(NRuns::from_bytes(&bytes, 1 << 30).is_err());
    }
}
