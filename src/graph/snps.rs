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
//! Editing SNPs towards the index before encoding.
//!
//! A contig that differs from the index by a SNP is encoded as several short
//! dictionary entries around the SNP. Replacing the contig's base with the
//! index's base before encoding gives fewer, longer entries; the original base
//! is stored separately and put back after decoding.
//!
//! The SNPs are called with [kbo::call] against the index. A SNP pays off if
//! few contigs share it: the graph already stores the nodes of a SNP that many
//! contigs share only once. SNPs are identified across contigs by the change
//! and [FLANK] bases on both sides.
//!
//! ## Stored SNPs
//!
//! A SNP that several contigs share is stored in full once. Each SNP gets an
//! id when it is stored in full ("new"), in colour order. A later contig with
//! the same SNP (same change and [FLANK] bases on both sides) stores it as a
//! *repeat*: its position is not stored, the decoder finds it by searching
//! the contig's edited sequence for the SNP's [CONTEXT]-base context taken
//! from the contig that has it as new. A SNP is only a repeat if its context
//! is the same as there and occurs once in the contig.
//!
//! Every colour has an id list: its new SNPs in position order, then its
//! repeats. A colour with repeats names a *reference*, an earlier colour, and
//! which entries of the reference's id list it shares; its repeats are those
//! entries. Related contigs share runs of SNPs, so this costs a few bits per
//! repeat. Repeats of SNPs not in the reference are stored as new.
//!
//! Per colour, in contig order, the bit stream has:
//!
//! - its new SNPs: the gap from the base after the previous new SNP (from the
//!   start of the contig for the first one) as a Huffman-coded bit length
//!   followed by the lower bits, and the original base as its distance from
//!   the index's base (1, 2 or 3 in A, C, G, T order; 2 is a transition), also
//!   Huffman coded. The index's base is read from the decoded (edited) contig.
//! - with a reference: one bit for how the shared entries are coded, then
//!   either a flag per entry of the reference's id list, or their number and
//!   the gaps between their indices, both Elias gamma coded.
//!
//! The number of new SNPs and the distance back to the reference (0 for
//! none) are stored per colour; the streams of all colours are compressed
//! together with LZMA.
//!
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use sbwt::LcsArray;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

use super::bits::{BitReader, BitWriter};
use super::huffman::HuffmanCode;
use super::lzma_wrapper::{self, LzmaData};

/// Default for the largest number of contigs that may share a SNP for it to
/// be edited.
pub const DEFAULT_MAX_SHARED: usize = 5;

/// Bases on both sides of a SNP that identify it across contigs.
pub const FLANK: usize = 32;

/// Length of the context that locates a repeated SNP: [FLANK] bases on both
/// sides of the index's base.
pub const CONTEXT: usize = 2 * FLANK + 1;

/// Number of gap symbols: the bit length of gap + 1, at most 64.
const N_GAP_SYMBOLS: usize = 65;

/// Number of base symbols: the distance from the index's base, 1..=3.
const N_BASE_SYMBOLS: usize = 4;

/// A SNP of a contig against the index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Snp {
    /// Position in the contig.
    pub pos: u32,

    /// The contig's base.
    pub base: u8,

    /// The index's base.
    pub index_base: u8,
}

/// A SNP and the hash that identifies it across contigs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CalledSnp {
    pub snp: Snp,
    pub key: u64,
}

fn base_code(
    base: u8,
) -> Option<u8> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Hash of the SNP with [FLANK] bases on both sides.
fn key(
    sequence: &[u8],
    snp: &Snp,
) -> u64 {
    let pos = snp.pos as usize;
    let start = pos.saturating_sub(FLANK);
    let end = (pos + 1 + FLANK).min(sequence.len());
    let mut hasher = DefaultHasher::new();
    sequence[start..pos].hash(&mut hasher);
    snp.base.hash(&mut hasher);
    snp.index_base.hash(&mut hasher);
    sequence[(pos + 1)..end].hash(&mut hasher);
    hasher.finish()
}

/// SNPs among the variants kbo calls for `variants`: sorted, without
/// overlapping calls, and with A, C, G or T on both sides.
fn snps_from_variants(
    sequence: &[u8],
    mut variants: Vec<kbo::variant_calling::Variant>,
) -> Vec<Snp> {
    variants.sort_by_key(|v| v.query_pos);
    let mut snps: Vec<Snp> = Vec::new();
    let mut prev_end = 0;
    let mut prev_pos = None;
    for v in variants {
        // Calls that overlap an earlier one are dropped, whatever their type.
        let end = v.query_pos + v.query_chars.len();
        if (v.query_chars.is_empty() && v.ref_chars.is_empty())
            || end > sequence.len()
            || v.query_pos < prev_end
            || prev_pos == Some(v.query_pos)
            || !v.ref_chars.iter().all(|&c| base_code(c).is_some()) {
            continue;
        }
        prev_end = end;
        prev_pos = Some(v.query_pos);
        if v.query_chars.len() == 1 && v.ref_chars.len() == 1 {
            let (base, index_base) = (sequence[v.query_pos], v.ref_chars[0]);
            if base_code(base).is_some() && base != index_base {
                snps.push(Snp { pos: v.query_pos.try_into().expect("contig position does not fit in u32"), base, index_base });
            }
        }
    }
    snps
}

/// Calls the SNPs of `sequence` against the index.
pub fn call_snps(
    sequence: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> Vec<CalledSnp> {
    let SbwtIndexVariant::SubsetMatrix(index) = sbwt;
    let mut opts = kbo::CallOpts::default();
    opts.sbwt_build_opts.k = index.k();
    opts.sbwt_build_opts.num_threads = 1;
    let variants = kbo::call(sbwt, lcs, sequence, opts);
    snps_from_variants(sequence, variants).into_iter().map(|snp| CalledSnp { snp, key: key(sequence, &snp) }).collect()
}

/// Error from reading the stored SNPs.
#[derive(Debug)]
pub struct SnpDecodeError(pub String);

impl std::fmt::Display for SnpDecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SNP edits: {}", self.0)
    }
}

impl std::error::Error for SnpDecodeError {}

/// An edited SNP with what is needed to store it as a repeat.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EditedSnp {
    pub snp: Snp,
    /// Identifies the SNP across contigs, see [FLANK].
    pub key: u64,
    /// The [CONTEXT] bases around the SNP in the edited contig (None at the
    /// contig ends).
    pub context: Option<[u8; CONTEXT]>,
    /// Whether the context occurs only once in the edited contig.
    pub unique: bool,
}

/// The [CONTEXT] bases of `sequence` centred on `pos` (None too close to an
/// end).
pub fn context_at(
    sequence: &[u8],
    pos: usize,
) -> Option<[u8; CONTEXT]> {
    if pos < FLANK || pos + FLANK + 1 > sequence.len() {
        return None
    }
    sequence[pos - FLANK..pos + FLANK + 1].try_into().ok()
}

const CONTEXT_HASH_BASE: u64 = 0x100000001B3;

fn context_hash(
    context: &[u8],
) -> u64 {
    context.iter().fold(0_u64, |h, &c| h.wrapping_mul(CONTEXT_HASH_BASE).wrapping_add(c as u64))
}

/// Start of every window of `sequence` that equals one of `contexts`, per
/// context (rolling hash, then a byte comparison).
fn find_contexts(
    sequence: &[u8],
    contexts: &[[u8; CONTEXT]],
) -> Vec<Vec<usize>> {
    let mut found: Vec<Vec<usize>> = vec![Vec::new(); contexts.len()];
    if contexts.is_empty() || sequence.len() < CONTEXT {
        return found
    }
    let mut by_hash: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, context) in contexts.iter().enumerate() {
        by_hash.entry(context_hash(context)).or_default().push(i);
    }
    let top = (0..CONTEXT - 1).fold(1_u64, |a, _| a.wrapping_mul(CONTEXT_HASH_BASE));
    let mut h = 0_u64;
    for i in 0..sequence.len() {
        if i >= CONTEXT {
            h = h.wrapping_sub((sequence[i - CONTEXT] as u64).wrapping_mul(top));
        }
        h = h.wrapping_mul(CONTEXT_HASH_BASE).wrapping_add(sequence[i] as u64);
        if i + 1 >= CONTEXT {
            if let Some(ids) = by_hash.get(&h) {
                let start = i + 1 - CONTEXT;
                for &id in ids {
                    if sequence[start..=i] == contexts[id] {
                        found[id].push(start);
                    }
                }
            }
        }
    }
    found
}

/// The SNPs of one contig with their contexts in `edited`, the contig with
/// the index's base at every SNP.
pub fn with_contexts(
    edited: &[u8],
    snps: Vec<CalledSnp>,
) -> Vec<EditedSnp> {
    let contexts: Vec<Option<[u8; CONTEXT]>> = snps.iter().map(|s| context_at(edited, s.snp.pos as usize)).collect();
    let present: Vec<[u8; CONTEXT]> = contexts.iter().flatten().copied().collect();
    let found = find_contexts(edited, &present);
    let mut counts = found.iter().map(|f| f.len());
    snps.into_iter().zip(contexts).map(|(s, context)| {
        let unique = context.is_some() && counts.next() == Some(1);
        EditedSnp { snp: s.snp, key: s.key, context, unique }
    }).collect()
}

/// Edited SNPs of one colour, as stored.
#[derive(serde::Serialize, serde::Deserialize, Debug, Default, Clone, PartialEq, Eq)]
pub struct SnpStream {
    /// Number of new SNPs in the stream.
    pub n_snps: u32,

    /// Colours back to the reference (0: no repeats).
    pub reference: u32,

    /// New SNPs and the shared entries of the reference, most significant
    /// bit first.
    pub bits: Vec<u8>,
}

/// Edited SNPs of every colour, see the module docs.
#[derive(serde::Serialize, serde::Deserialize, Debug, Default, Clone, PartialEq, Eq)]
struct StoredSnps {
    /// Huffman code lengths of the gap symbols.
    gap_code: Vec<u8>,

    /// Huffman code lengths of the base symbols.
    base_code: Vec<u8>,

    /// One stream per colour.
    streams: Vec<SnpStream>,
}

/// SNPs of one colour, decoded.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
struct ColourSnps {
    /// New SNPs: position and the base's distance from the index's base.
    new: Vec<(u32, u8)>,
    /// Ids of the repeated SNPs.
    repeats: Vec<u32>,
}

/// Edited SNPs of every colour.
///
/// Serialized as LZMA-compressed [StoredSnps].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnpEdits {
    gap_code: HuffmanCode,
    base_code: HuffmanCode,
    streams: Vec<SnpStream>,
    colours: Vec<ColourSnps>,
    /// Colour and index in its new SNPs of every id.
    ids: Vec<(u32, u32)>,
}

fn gap_symbol(
    gap: u64,
) -> usize {
    (64 - (gap + 1).leading_zeros()) as usize
}

/// Distance of `base` from `index_base`, 1..=3.
fn base_symbol(
    base: u8,
    index_base: u8,
) -> usize {
    ((base_code(base).unwrap() + 4 - base_code(index_base).unwrap()) % 4) as usize
}

fn write_gamma(
    bits: &mut BitWriter,
    value: u64,
) {
    // value >= 1
    let n = 64 - value.leading_zeros();
    bits.write_bits(0, n - 1);
    bits.write_bits(value, n);
}

fn read_gamma(
    bits: &mut BitReader,
) -> Option<u64> {
    let mut zeros = 0;
    while bits.read_bits(1)? == 0 {
        zeros += 1;
        if zeros > 63 { return None }
    }
    Some((1_u64 << zeros) | bits.read_bits(zeros)?)
}

impl SnpEdits {
    /// Codes the SNPs of every colour, each sorted by position, all as new.
    pub fn new(
        snps: &[Vec<Snp>],
    ) -> Self {
        let edited: Vec<Vec<EditedSnp>> = snps.iter().map(|colour| colour.iter().map(|&snp| EditedSnp { snp, key: 0, context: None, unique: false }).collect()).collect();
        SnpEdits::with_repeats(&edited)
    }

    /// Codes the SNPs of every colour, each sorted by position, storing
    /// shared SNPs once (see the module docs).
    pub fn with_repeats(
        snps: &[Vec<EditedSnp>],
    ) -> Self {
        // First occurrence (id) and context of every SNP stored as new.
        let mut first: HashMap<u64, (u32, [u8; CONTEXT])> = HashMap::new();
        // Colours whose id list holds an id.
        let mut holders: HashMap<u32, Vec<u32>> = HashMap::new();
        let mut lists: Vec<Vec<u32>> = Vec::with_capacity(snps.len());
        let mut next_id = 0_u32;
        // Per colour: new SNPs, reference distance and the shared entries' indices.
        let mut plan: Vec<(Vec<Snp>, u32, Vec<u32>)> = Vec::with_capacity(snps.len());
        for (c, colour) in snps.iter().enumerate() {
            let c = c as u32;
            let candidate: Vec<Option<u32>> = colour.iter().map(|s| {
                if !s.unique { return None }
                let (id, context) = first.get(&s.key)?;
                (s.context.as_ref() == Some(context)).then_some(*id)
            }).collect();
            let mut votes: HashMap<u32, usize> = HashMap::new();
            for id in candidate.iter().flatten() {
                for &r in holders.get(id).map_or(&[][..], |v| v.as_slice()) { *votes.entry(r).or_default() += 1; }
            }
            let reference = votes.into_iter().max_by_key(|&(r, v)| (v, r)).map(|(r, _)| r);
            let mut shared: Vec<u32> = Vec::new();
            let mut repeat_ids: std::collections::HashSet<u32> = std::collections::HashSet::new();
            if let Some(r) = reference {
                let wanted: std::collections::HashSet<u32> = candidate.iter().flatten().copied().collect();
                for (i, id) in lists[r as usize].iter().enumerate() {
                    if wanted.contains(id) { shared.push(i as u32); repeat_ids.insert(*id); }
                }
            }
            let mut list: Vec<u32> = Vec::new();
            let mut new: Vec<Snp> = Vec::new();
            for (s, cand) in colour.iter().zip(&candidate) {
                if cand.is_some_and(|id| repeat_ids.contains(&id)) { continue }
                new.push(s.snp);
                if let (Some(context), false) = (s.context, first.contains_key(&s.key)) {
                    first.insert(s.key, (next_id, context));
                }
                list.push(next_id);
                next_id += 1;
            }
            if let Some(r) = reference {
                list.extend(shared.iter().map(|&i| lists[r as usize][i as usize]));
            }
            for &id in &list { holders.entry(id).or_default().push(c); }
            lists.push(list);
            let distance = if shared.is_empty() { 0 } else { c - reference.unwrap() };
            plan.push((new, distance, if distance == 0 { Vec::new() } else { shared }));
        }

        let mut gap_counts = vec![0_u64; N_GAP_SYMBOLS];
        let mut base_counts = vec![0_u64; N_BASE_SYMBOLS];
        for (new, _, _) in &plan {
            let mut next = 0_u64;
            for snp in new {
                gap_counts[gap_symbol(snp.pos as u64 - next)] += 1;
                base_counts[base_symbol(snp.base, snp.index_base)] += 1;
                next = snp.pos as u64 + 1;
            }
        }
        let gap_code = HuffmanCode::from_counts(&gap_counts);
        let base_code = HuffmanCode::from_counts(&base_counts);

        let streams: Vec<SnpStream> = plan.iter().enumerate().map(|(c, (new, distance, shared))| {
            let mut bits = BitWriter::new();
            let mut next = 0_u64;
            for snp in new {
                assert!(snp.pos as u64 >= next, "SNPs are not sorted by position");
                let gap = snp.pos as u64 - next;
                let symbol = gap_symbol(gap);
                gap_code.write(symbol, &mut bits);
                bits.write_bits(gap + 1, symbol as u32 - 1);
                base_code.write(base_symbol(snp.base, snp.index_base), &mut bits);
                next = snp.pos as u64 + 1;
            }
            if *distance > 0 {
                let n_ref = lists[c - *distance as usize].len();
                // a flag per entry, or the count and the gaps (Elias gamma)
                let gaps_len: u64 = {
                    let mut prev = 0_u64;
                    let mut len = 2 * (64 - (shared.len() as u64).leading_zeros() as u64) - 1;
                    for &i in shared.iter() { let g = i as u64 - prev + 1; len += 2 * (64 - g.leading_zeros() as u64) - 1; prev = i as u64 + 1; }
                    len
                };
                if (n_ref as u64) <= gaps_len {
                    bits.write_bits(0, 1);
                    let mut j = 0;
                    for i in 0..n_ref as u32 {
                        let flag = j < shared.len() && shared[j] == i;
                        if flag { j += 1; }
                        bits.write_bits(flag as u64, 1);
                    }
                } else {
                    bits.write_bits(1, 1);
                    write_gamma(&mut bits, shared.len() as u64);
                    let mut prev = 0_u64;
                    for &i in shared.iter() { write_gamma(&mut bits, i as u64 - prev + 1); prev = i as u64 + 1; }
                }
            }
            SnpStream { n_snps: new.len().try_into().expect("too many SNPs"), reference: *distance, bits: bits.finish() }
        }).collect();

        SnpEdits::from_streams(gap_code, base_code, streams).expect("SNP streams written above decode")
    }

    /// Decodes the streams into new SNPs and repeats per colour.
    fn from_streams(
        gap_code: HuffmanCode,
        base_code: HuffmanCode,
        streams: Vec<SnpStream>,
    ) -> Result<Self, SnpDecodeError> {
        let mut colours: Vec<ColourSnps> = Vec::with_capacity(streams.len());
        let mut lists: Vec<Vec<u32>> = Vec::with_capacity(streams.len());
        let mut ids: Vec<(u32, u32)> = Vec::new();
        for (c, stream) in streams.iter().enumerate() {
            let truncated = || SnpDecodeError(format!("SNP stream of colour {} is truncated", c));
            let mut bits = BitReader::new(&stream.bits);
            let mut new: Vec<(u32, u8)> = Vec::with_capacity(stream.n_snps as usize);
            let mut list: Vec<u32> = Vec::new();
            let mut next = 0_u64;
            for i in 0..stream.n_snps {
                let symbol = gap_code.read(&mut bits).ok_or_else(truncated)?;
                if symbol == 0 {
                    return Err(truncated())
                }
                let gap = ((1_u64 << (symbol - 1)) | bits.read_bits(symbol as u32 - 1).ok_or_else(truncated)?) - 1;
                let pos = next + gap;
                let base = base_code.read(&mut bits).ok_or_else(truncated)?;
                if base == 0 {
                    return Err(SnpDecodeError(format!("SNP of colour {} at {} does not change the base", c, pos)))
                }
                let pos: u32 = pos.try_into().map_err(|_| truncated())?;
                new.push((pos, base as u8));
                list.push(ids.len() as u32);
                ids.push((c as u32, i));
                next = pos as u64 + 1;
            }
            let mut repeats: Vec<u32> = Vec::new();
            if stream.reference > 0 {
                let r = c.checked_sub(stream.reference as usize).ok_or_else(|| SnpDecodeError(format!("SNP reference of colour {} is before the first colour", c)))?;
                let reference = &lists[r];
                if bits.read_bits(1).ok_or_else(truncated)? == 0 {
                    for &id in reference.iter() {
                        if bits.read_bits(1).ok_or_else(truncated)? == 1 { repeats.push(id); }
                    }
                } else {
                    let n = read_gamma(&mut bits).ok_or_else(truncated)?;
                    let mut prev = 0_u64;
                    for _ in 0..n {
                        let i = prev + read_gamma(&mut bits).ok_or_else(truncated)? - 1;
                        repeats.push(*reference.get(i as usize).ok_or_else(truncated)?);
                        prev = i + 1;
                    }
                }
            }
            list.extend_from_slice(&repeats);
            lists.push(list);
            colours.push(ColourSnps { new, repeats });
        }
        Ok(SnpEdits { gap_code, base_code, streams, colours, ids })
    }

    /// Number of colours.
    pub fn n_colors(&self) -> usize {
        self.colours.len()
    }

    /// Total number of edited SNPs, new and repeated.
    pub fn n_snps(&self) -> usize {
        self.colours.iter().map(|c| c.new.len() + c.repeats.len()).sum()
    }

    /// Number of SNPs stored as repeats.
    pub fn n_repeats(&self) -> usize {
        self.colours.iter().map(|c| c.repeats.len()).sum()
    }

    /// Colours and positions whose contexts (see [context_at], in the edited
    /// contig of that colour) [SnpEdits::revert] needs for `color`.
    pub fn context_sources(
        &self,
        color: u32,
    ) -> Vec<(u32, u32)> {
        self.colours.get(color as usize).map_or_else(Vec::new, |c| c.repeats.iter().map(|&id| {
            let (f, i) = self.ids[id as usize];
            (f, self.colours[f as usize].new[i as usize].0)
        }).collect())
    }

    /// Puts the original bases of `color` back into its decoded contig.
    /// `context` gives the context at a position of an earlier colour's edited
    /// contig (for the pairs in [SnpEdits::context_sources]).
    pub fn revert(
        &self,
        color: u32,
        sequence: &mut [u8],
        context: &dyn Fn(u32, u32) -> Option<[u8; CONTEXT]>,
    ) -> Result<(), SnpDecodeError> {
        let snps = self.colours.get(color as usize).ok_or_else(|| SnpDecodeError(format!("no SNP stream for colour {}", color)))?;
        // Repeats are located in the edited sequence before any base changes.
        let mut targets: Vec<[u8; CONTEXT]> = Vec::with_capacity(snps.repeats.len());
        let mut distances: Vec<u8> = Vec::with_capacity(snps.repeats.len());
        for &id in &snps.repeats {
            let (f, i) = self.ids[id as usize];
            let (pos, distance) = self.colours[f as usize].new[i as usize];
            targets.push(context(f, pos).ok_or_else(|| SnpDecodeError(format!("no context at {} of colour {} for colour {}", pos, f, color)))?);
            distances.push(distance);
        }
        let found = find_contexts(sequence, &targets);
        let mut changes: Vec<(usize, u8)> = Vec::with_capacity(snps.new.len() + targets.len());
        for (j, starts) in found.iter().enumerate() {
            if starts.len() != 1 {
                return Err(SnpDecodeError(format!("repeated SNP {} of colour {} found {} times", j, color, starts.len())))
            }
            changes.push((starts[0] + FLANK, distances[j]));
        }
        changes.extend(snps.new.iter().map(|&(pos, d)| (pos as usize, d)));
        for (pos, distance) in changes {
            let index_base = sequence.get(pos).and_then(|&b| base_code(b))
                .ok_or_else(|| SnpDecodeError(format!("SNP of colour {} at {} is not on an A, C, G or T of the contig", color, pos)))?;
            sequence[pos] = b"ACGT"[(index_base as usize + distance as usize) % 4];
        }
        Ok(())
    }

    fn to_stored(
        &self,
    ) -> Result<LzmaData, String> {
        let stored = StoredSnps {
            gap_code: self.gap_code.lengths().to_vec(),
            base_code: self.base_code.lengths().to_vec(),
            streams: self.streams.clone(),
        };
        let bytes = postcard::to_allocvec(&stored).map_err(|e| format!("SNPs: {}", e))?;
        lzma_wrapper::compress(&bytes).map_err(|e| format!("SNPs: {}", e))
    }

    fn from_stored(
        stored: LzmaData,
    ) -> Result<Self, String> {
        let bytes = lzma_wrapper::decompress(&stored).map_err(|e| format!("SNPs: {}", e))?;
        let stored: StoredSnps = postcard::from_bytes(&bytes).map_err(|e| format!("SNPs: {}", e))?;
        if stored.gap_code.len() != N_GAP_SYMBOLS || stored.base_code.len() != N_BASE_SYMBOLS {
            return Err("SNPs: invalid code tables".to_string())
        }
        SnpEdits::from_streams(
            HuffmanCode::from_lengths(stored.gap_code).map_err(|e| format!("SNP gaps: {}", e))?,
            HuffmanCode::from_lengths(stored.base_code).map_err(|e| format!("SNP bases: {}", e))?,
            stored.streams,
        ).map_err(|e| e.to_string())
    }
}

impl serde::Serialize for SnpEdits {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.to_stored().map_err(serde::ser::Error::custom)?.serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for SnpEdits {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let stored = LzmaData::deserialize(deserializer)?;
        SnpEdits::from_stored(stored).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kbo::variant_calling::Variant;

    fn snp(pos: u32, base: u8, index_base: u8) -> Snp {
        Snp { pos, base, index_base }
    }

    /// `sequence` with the index's base at every SNP.
    fn edit(
        sequence: &[u8],
        snps: &[Snp],
    ) -> Vec<u8> {
        let mut edited = sequence.to_vec();
        for snp in snps {
            edited[snp.pos as usize] = snp.index_base;
        }
        edited
    }

    #[test]
    fn edit_and_revert_roundtrip() {
        let original = b"ACGTACGTACGTTTGCA".to_vec();
        let colours = vec![
            vec![snp(0, b'A', b'G'), snp(1, b'C', b'T'), snp(7, b'T', b'A'), snp(16, b'A', b'C')],
            vec![],
            vec![snp(5, b'C', b'G')],
        ];
        let edits: SnpEdits = postcard::from_bytes(&postcard::to_allocvec(&SnpEdits::new(&colours)).unwrap()).unwrap();
        assert_eq!(edits.n_colors(), 3);
        assert_eq!(edits.n_snps(), 5);
        for (color, snps) in colours.iter().enumerate() {
            let mut sequence = edit(&original, snps);
            if !snps.is_empty() {
                assert_ne!(sequence, original);
            }
            edits.revert(color as u32, &mut sequence, &|_, _| None).unwrap();
            assert_eq!(sequence, original);
        }
        assert!(edits.revert(3, &mut original.clone(), &|_, _| None).is_err());
    }

    #[test]
    fn corrupt_streams_are_errors() {
        let edits = SnpEdits::new(&[vec![snp(3, b'A', b'C'), snp(1000, b'G', b'A')]]);
        // The contig is too short for the second SNP.
        assert!(edits.revert(0, &mut b"ACGTACGT".to_vec(), &|_, _| None).is_err());
        let mut streams = edits.streams.clone();
        streams[0].bits.truncate(1);
        assert!(SnpEdits::from_streams(edits.gap_code.clone(), edits.base_code.clone(), streams).is_err());
    }

    #[test]
    fn shared_snps_are_stored_once() {
        // Three contigs share two SNPs; the third also has one of its own. The
        // later two store the shared SNPs as repeats of the first contig's.
        let original = random_sequence(2000, 11);
        let shared = [snp(300, original[300], b"ACGT"[(base_code(original[300]).unwrap() as usize + 1) % 4]),
                      snp(900, original[900], b"ACGT"[(base_code(original[900]).unwrap() as usize + 2) % 4])];
        let own = snp(1500, original[1500], b"ACGT"[(base_code(original[1500]).unwrap() as usize + 3) % 4]);
        let colours: Vec<Vec<Snp>> = vec![shared.to_vec(), shared.to_vec(), vec![shared[0], shared[1], own]];
        let contigs: Vec<Vec<u8>> = colours.iter().map(|c| original.clone()).collect();
        let edited_seqs: Vec<Vec<u8>> = colours.iter().zip(&contigs).map(|(c, s)| edit(s, c)).collect();
        let called: Vec<Vec<EditedSnp>> = colours.iter().zip(&contigs).zip(&edited_seqs).map(|((c, s), e)| {
            with_contexts(e, c.iter().map(|&snp| CalledSnp { snp, key: key(s, &snp) }).collect())
        }).collect();
        let edits: SnpEdits = postcard::from_bytes(&postcard::to_allocvec(&SnpEdits::with_repeats(&called)).unwrap()).unwrap();
        assert_eq!(edits.n_snps(), 7);
        assert_eq!(edits.n_repeats(), 4);
        assert!(edits.context_sources(0).is_empty());
        assert_eq!(edits.context_sources(2), vec![(0, 300), (0, 900)]);
        let context = |f: u32, pos: u32| context_at(&edited_seqs[f as usize], pos as usize);
        for (color, contig) in contigs.iter().enumerate() {
            let mut sequence = edited_seqs[color].clone();
            edits.revert(color as u32, &mut sequence, &context).unwrap();
            assert_eq!(&sequence, contig, "colour {}", color);
        }
        // Without the contexts the repeats cannot be put back.
        assert!(edits.revert(1, &mut edited_seqs[1].clone(), &|_, _| None).is_err());
    }

    fn random_sequence(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len).map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize]
        }).collect()
    }

    #[test]
    fn builder_edits_snps_and_decoding_puts_them_back() {
        use super::super::{GraphBuilder, SuffixLen, Walkable, decode_path};

        let reference = random_sequence(4000, 7);
        // kbo calls SNPs only if k is well above the length of random matches.
        let k = 63;
        let mut opts = kbo::BuildOpts::default();
        opts.k = k as usize;
        opts.add_revcomp = true;
        opts.build_select = true;
        let (sbwt, lcs) = kbo::build(&[reference.clone()], opts);

        // Two contigs with the same SNPs and one with its own, more than k apart.
        let with_snps = |positions: &[usize]| {
            let mut contig = reference.clone();
            for &i in positions {
                contig[i] = b"ACGT"[(base_code(contig[i]).unwrap() as usize + 1) % 4];
            }
            contig
        };
        let contigs = [with_snps(&[500, 1500, 2500]), with_snps(&[500, 1500, 2500]), with_snps(&[1000, 3000])];

        for max_shared in [1, 2] {
            let mut builder = GraphBuilder::new().collapse_chains(true).edit_snps(true).edit_snps_max_shared(max_shared);
            assert!(builder.add_contig(&contigs[0], &sbwt, &lcs).is_err(), "SNPs must be called first");
            for contig in &contigs {
                builder.call_snps(contig, &sbwt, &lcs);
            }
            for contig in &contigs {
                builder.add_contig(contig, &sbwt, &lcs).unwrap();
            }
            let edits = builder.snp_edits().unwrap();
            // The shared SNPs are only edited if two contigs may share a SNP.
            assert_eq!(edits.n_snps(), if max_shared == 1 { 2 } else { 8 });
            let (graph, choices) = builder.finish();
            // Decode every edited contig first: repeated SNPs need the contexts of
            // earlier contigs.
            let mut edited_seqs: Vec<Vec<u8>> = Vec::new();
            for color in 0..contigs.len() {
                let (len, dictionary) = graph.walk_path(&choices, color as u32).unwrap();
                let dictionary = dictionary.into_iter().map(|(colex, len)| (colex, if len == 0 {
                    SuffixLen::Implicit
                } else if len > k {
                    SuffixLen::Explicit { suffix_len: len }
                } else {
                    SuffixLen::Plain { suffix_len: len }
                })).collect();
                edited_seqs.push(decode_path(dictionary, len as usize, &sbwt, Some(&lcs)).unwrap());
            }
            if max_shared > 1 {
                assert!(edits.n_repeats() > 0, "the second contig repeats the first one's SNPs");
            }
            for (color, contig) in contigs.iter().enumerate() {
                let mut decoded = edited_seqs[color].clone();
                edits.revert(color as u32, &mut decoded, &|f, pos| context_at(&edited_seqs[f as usize], pos as usize)).unwrap();
                assert_eq!(&decoded, contig, "colour {}", color);
            }
        }
    }

    #[test]
    fn only_non_overlapping_snps_are_kept() {
        let sequence = b"AAAAACCCCCGGGGGTTTTT".to_vec();
        let variant = |pos: usize, query: &[u8], reference: &[u8]| Variant { query_pos: pos, query_chars: query.to_vec(), ref_chars: reference.to_vec() };
        let variants = vec![
            variant(12, b"G", b"T"),   // SNP
            variant(2, b"AAA", b""),   // deletion from the index, covers 2..5
            variant(3, b"A", b"C"),    // overlaps the deletion
            variant(12, b"G", b"A"),   // second call at the same position
            variant(16, b"T", b"N"),   // not a base in the index
            variant(7, b"C", b"GG"),   // not a SNP
            variant(9, b"C", b"A"),    // SNP
        ];
        assert_eq!(snps_from_variants(&sequence, variants), vec![snp(9, b'C', b'A'), snp(12, b'G', b'T')]);
    }
}
