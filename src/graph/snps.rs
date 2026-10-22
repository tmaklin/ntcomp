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
//! Editing SNPs and small indels towards the index before encoding.
//!
//! A contig that differs from the index by a SNP is encoded as several short
//! dictionary entries around the SNP. Replacing the contig's base with the
//! index's base before encoding gives fewer, longer entries; the original base
//! is stored separately and put back after decoding. Small indels are edited
//! the same way: contig bases the index does not have are removed (a
//! deletion, [Change::Del]) and index bases the contig does not have are
//! added (an insertion, [Change::Ins]). Below, "SNP" means any of these edits.
//!
//! The edits are called by walking the index, see [call_snps].
//!
//! ## Stored SNPs
//!
//! A SNP that several contigs share is stored in full once. Each SNP gets an
//! id when it is stored in full ("new"), in colour order. A later contig with
//! the same SNP (the same change and [CONTEXT]-base context in the edited
//! contig, centred on the SNP's position there) stores it as a *repeat*: its
//! position is not stored, the decoder finds it by searching the contig's
//! edited sequence for the context taken from the contig that has it as new.
//! A SNP is only a repeat if its context occurs once in the contig.
//!
//! Every colour has an id list: its new SNPs in position order, then its
//! repeats. A colour with repeats may name a *reference*, the earlier colour
//! holding most of them, and which entries of the reference's id list it
//! shares; related contigs share runs of SNPs, so this costs a fraction of a
//! bit per entry. Its other repeats are *direct*: their ids.
//!
//! All colours are coded as one stream with the arithmetic coder
//! ([crate::coder::arith]): the number of colours, then per colour
//! - the number of new SNPs (adaptive Elias gamma, by whether the previous
//!   colour had none);
//! - its new SNPs: the gap from the base after the previous new SNP (from the
//!   start of the contig for the first one), adaptive Elias gamma by the
//!   bit length of the previous gap, and its change as a symbol of an
//!   adaptive 2-bit tree: 0–2 a substitution, the original base's distance
//!   from the index's base minus one (1, 2 or 3 in A, C, G, T order; 2 is a
//!   transition; the index's base is read from the decoded, edited, contig),
//!   3 an indel: an adaptive flag (deletion or insertion), the length minus
//!   one (an adaptive 2-bit tree per kind) and for a deletion its bases (an
//!   adaptive 2-bit tree). The position is in the edited contig: where a
//!   deletion's bases are put back, the first of an insertion's bases;
//! - whether it has a reference; if so the distance back to it and a flag per
//!   entry of the reference's id list (adaptive, by the previous two flags,
//!   whether the entry is one of the reference's new SNPs and whether it is
//!   the last);
//! - the number of direct repeats and their ids in increasing order: the
//!   first as the distance back from the newest id of the earlier colours,
//!   the others as differences (adaptive Elias gamma, by the bit length of
//!   the previous one).
//!
use std::collections::{HashMap, HashSet};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use sbwt::{ExtendRight, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};
use sbwt::sbwt_index_variant::SbwtIndexVariant;

use crate::coder::arith::{code_bit, code_tree, Coder, Decoder, Encoder, GammaModel, C0};

/// Default for the largest number of contigs that may share a SNP for it to
/// be edited: no limit, every SNP is edited.
pub const DEFAULT_MAX_SHARED: usize = usize::MAX;

/// Bases on both sides of a SNP that identify it across contigs.
pub const FLANK: usize = 32;

/// Length of the context that locates a repeated SNP: [FLANK] bases on both
/// sides of the index's base.
pub const CONTEXT: usize = 2 * FLANK + 1;

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

/// A change of a contig against the index.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Change {
    /// A SNP: the contig's base `base` where the index has `index_base`.
    Sub { base: u8, index_base: u8 },
    /// Contig bases the index does not have, removed when editing.
    Del { bases: Vec<u8> },
    /// Index bases the contig does not have, added when editing.
    Ins { bases: Vec<u8> },
}

/// A change at a position of a contig: the first base it changes or removes,
/// or the base an insertion goes before.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Edit {
    pub pos: u32,
    pub change: Change,
}

impl From<Snp> for Edit {
    fn from(snp: Snp) -> Self {
        Edit { pos: snp.pos, change: Change::Sub { base: snp.base, index_base: snp.index_base } }
    }
}

/// A called edit (in the contig's coordinates) and the hash that identifies
/// it across contigs ([FLANK] contig bases on both sides and the change).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CalledSnp {
    pub edit: Edit,
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

/// Hash of the edit with [FLANK] bases of `sequence` on both sides (before
/// its position and after the bases it changes or removes).
fn key(
    sequence: &[u8],
    edit: &Edit,
) -> u64 {
    let pos = edit.pos as usize;
    let after = match &edit.change {
        Change::Sub { .. } => pos + 1,
        Change::Del { bases } => pos + bases.len(),
        Change::Ins { .. } => pos,
    };
    let mut hasher = DefaultHasher::new();
    sequence[pos.saturating_sub(FLANK)..pos].hash(&mut hasher);
    edit.change.hash(&mut hasher);
    sequence[after..(after + FLANK).min(sequence.len())].hash(&mut hasher);
    hasher.finish()
}

/// Most edits in one stretch of the contig not in the index, see
/// [call_snps].
const MAX_EDITS: usize = 5;

/// Most contig bases one deletion removes, see [call_snps].
const MAX_DELETED: usize = 3;

/// Most index bases one insertion adds, see [call_snps].
const MAX_INSERTED: usize = 2;

/// A pattern in the index: its colex interval and length (at most _k_).
type Pattern = (std::ops::Range<usize>, usize);

/// One step right in the index: the interval of `pattern` followed by `c`,
/// dropping the pattern's first character first if it is _k_ long. None if
/// that is not in the index.
fn step_right(
    sbwt: &SbwtIndex<SubsetMatrix>,
    lcs: &LcsArray,
    pattern: &Pattern,
    c: u8,
) -> Option<Pattern> {
    let k = sbwt.k();
    let (mut interval, mut len) = pattern.clone();
    if len == k {
        // The k-mers sharing the last k - 1 characters.
        while interval.start > 0 && lcs.access(interval.start) >= k - 1 {
            interval.start -= 1;
        }
        while interval.end < sbwt.n_sets() && lcs.access(interval.end) >= k - 1 {
            interval.end += 1;
        }
        len = k - 1;
    }
    let next = sbwt.extend_right(interval, c);
    (!next.is_empty()).then_some((next, len + 1))
}

/// Steps right in the index by the characters of `sequence` from `from`,
/// stopping before `until` or at the first one not in the index. Returns the
/// pattern and the first position not taken.
fn walk_right(
    sbwt: &SbwtIndex<SubsetMatrix>,
    lcs: &LcsArray,
    sequence: &[u8],
    mut pattern: Pattern,
    from: usize,
    until: usize,
) -> (Pattern, usize) {
    let mut pos = from;
    while pos < until {
        match step_right(sbwt, lcs, &pattern, sequence[pos]) {
            Some(next) => {
                pattern = next;
                pos += 1;
            },
            None => break,
        }
    }
    (pattern, pos)
}

/// An edit tried in a stretch the contig cannot extend: its change (as
/// the base, the number of contig bases removed, or the bases added), the
/// pattern after the walk following it, the contig position the walk
/// reached and whether that completes the stretch.
type Tried = (TriedChange, Pattern, usize, bool);

#[derive(Clone)]
enum TriedChange {
    Sub(u8),
    Del(usize),
    Ins(Vec<u8>),
}

/// Edits from `pos` until a walk takes _k_ contig bases from the last one,
/// starting from the match `pattern`: at most [MAX_EDITS], substitutions only
/// or also indels. None if the walk does not get through.
fn repair(
    sbwt: &SbwtIndex<SubsetMatrix>,
    lcs: &LcsArray,
    sequence: &[u8],
    mut pattern: Pattern,
    mut pos: usize,
    with_indels: bool,
) -> Option<(Vec<(usize, TriedChange)>, Pattern, usize)> {
    let (k, n) = (sbwt.k(), sequence.len());
    let mut found: Vec<(usize, TriedChange)> = Vec::new();
    loop {
        if found.len() == MAX_EDITS || pos + k > n || base_code(sequence[pos]).is_none() {
            return None
        }
        let mut tried: Vec<Tried> = Vec::new();
        for &b in b"ACGT" {
            if b == sequence[pos] {
                continue
            }
            let Some(after) = step_right(sbwt, lcs, &pattern, b) else { continue };
            let (walked, reached) = walk_right(sbwt, lcs, sequence, after, pos + 1, pos + k);
            tried.push((TriedChange::Sub(b), walked, reached, reached == pos + k));
        }
        if with_indels && !tried.iter().any(|t| t.3) {
            // Contig bases removed: the walk goes on after them.
            for d in 1..=MAX_DELETED {
                if pos + d + k > n || !sequence[pos..pos + d].iter().all(|&c| base_code(c).is_some()) {
                    break
                }
                let (walked, reached) = walk_right(sbwt, lcs, sequence, pattern.clone(), pos + d, pos + d + k);
                if reached > pos + d {
                    tried.push((TriedChange::Del(d), walked, reached, reached == pos + d + k));
                }
            }
            // Index bases added: every path of up to MAX_INSERTED steps.
            let mut paths: Vec<(Vec<u8>, Pattern)> = vec![(Vec::new(), pattern.clone())];
            for _ in 0..MAX_INSERTED {
                let mut longer = Vec::new();
                for (bases, path) in &paths {
                    for &b in b"ACGT" {
                        let Some(after) = step_right(sbwt, lcs, path, b) else { continue };
                        let mut bases = bases.clone();
                        bases.push(b);
                        let (walked, reached) = walk_right(sbwt, lcs, sequence, after.clone(), pos, pos + k);
                        if reached > pos {
                            tried.push((TriedChange::Ins(bases.clone()), walked, reached, reached == pos + k));
                        }
                        longer.push((bases, after));
                    }
                }
                paths = longer;
            }
        }
        // The first edit that completes the stretch (substitutions first,
        // then deletions and insertions by length), else the one that walks
        // furthest.
        let best = tried.iter().position(|t| t.3).or_else(|| {
            let furthest = tried.iter().map(|t| t.2).max()?;
            tried.iter().position(|t| t.2 == furthest)
        })?;
        let (change, walked, reached, complete) = tried.swap_remove(best);
        found.push((pos, change));
        if complete {
            return Some((found, walked, reached))
        }
        pattern = walked;
        pos = reached;
    }
}

/// Calls the edits (SNPs and small indels) of `sequence` against the index
/// by walking its de Bruijn graph; sorted by position.
///
/// Where the contig's match in the index (its matching statistics), at least
/// as long as a random match is unlikely to be ([random_match_threshold](kbo::derandomize::random_match_threshold)),
/// cannot be extended by the next base at p, every other base b is tried:
/// the match followed by b and then by the contig's bases, one step right in
/// the index at a time. If that walk takes the next _k_ − 1 bases, every
/// _k_-mer through p is in the index once p is b: a SNP. If no b gets that
/// far, the one that walks furthest is taken and the base where its walk
/// stops is substituted the same way, at most [MAX_EDITS] times. If
/// substitutions do not get the walk through, the same search is tried with
/// indels too: removing up to [MAX_DELETED] contig bases or adding up to
/// [MAX_INSERTED] index bases (every path of the index after the match);
/// a stretch whose walk does not get through has no edits. After an edit the
/// edited contig is walked on until its matching statistics are full _k_-mers
/// again, and a base not in the index there is repaired from that walk.
pub fn call_snps(
    sequence: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> Vec<CalledSnp> {
    let SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt;
    let k = sbwt.k();
    let n = sequence.len();
    if n <= k {
        return Vec::new()
    }
    let min_match = kbo::derandomize::random_match_threshold(k, sbwt.n_kmers(), 4, kbo::CallOpts::default().max_error_prob);
    let ms = StreamingIndex::new(sbwt, lcs).matching_statistics(sequence);

    let mut edits: Vec<Edit> = Vec::new();
    // The edited contig's match ending at i - 1, if a walk after an edit
    // stopped at i.
    let mut carried: Option<Pattern> = None;
    let mut i = 1;
    while i < n {
        let pattern = match carried.take() {
            Some(pattern) => pattern,
            None => {
                let (prev_len, len) = (ms[i - 1].0, ms[i].0);
                if prev_len < min_match || len >= (prev_len + 1).min(k) || base_code(sequence[i]).is_none() {
                    i += 1;
                    continue
                }
                (ms[i - 1].1.clone(), prev_len)
            },
        };

        let repaired = repair(sbwt, lcs, sequence, pattern.clone(), i, false)
            .or_else(|| repair(sbwt, lcs, sequence, pattern, i, true));
        match repaired {
            Some((found, mut pattern, mut pos)) => {
                for (at, change) in found {
                    let pos = at.try_into().expect("contig position does not fit in u32");
                    edits.push(match change {
                        TriedChange::Sub(b) => Edit { pos, change: Change::Sub { base: sequence[at], index_base: b } },
                        TriedChange::Del(d) => Edit { pos, change: Change::Del { bases: sequence[at..at + d].to_vec() } },
                        TriedChange::Ins(bases) => Edit { pos, change: Change::Ins { bases } },
                    });
                }
                // The original matching statistics still contain the edits here.
                while pos < n && ms[pos].0 < k {
                    match step_right(sbwt, lcs, &pattern, sequence[pos]) {
                        Some(next) => {
                            pattern = next;
                            pos += 1;
                        },
                        None => break,
                    }
                }
                if pos < n && ms[pos].0 < k {
                    carried = Some(pattern);
                }
                i = pos;
            },
            None => i += 1,
        }
    }
    edits.into_iter().map(|edit| CalledSnp { key: key(sequence, &edit), edit }).collect()
}

/// `sequence` with `edits` (in its coordinates, increasing and not
/// overlapping) applied, and the edits at their positions in the result.
/// None if an edit does not match the sequence.
pub fn apply_edits(
    sequence: &[u8],
    edits: Vec<Edit>,
) -> Option<(Vec<u8>, Vec<Edit>)> {
    let mut out: Vec<u8> = Vec::with_capacity(sequence.len() + 16);
    let mut moved: Vec<Edit> = Vec::with_capacity(edits.len());
    let mut from = 0;
    for edit in edits {
        let pos = edit.pos as usize;
        if pos < from || pos > sequence.len() {
            return None
        }
        out.extend_from_slice(&sequence[from..pos]);
        let at = out.len() as u32;
        match &edit.change {
            Change::Sub { base, index_base } => {
                if sequence.get(pos) != Some(base) {
                    return None
                }
                out.push(*index_base);
                from = pos + 1;
            },
            Change::Del { bases } => {
                if sequence.get(pos..pos + bases.len()) != Some(bases.as_slice()) {
                    return None
                }
                from = pos + bases.len();
            },
            Change::Ins { bases } => {
                out.extend_from_slice(bases);
                from = pos;
            },
        }
        moved.push(Edit { pos: at, change: edit.change });
    }
    out.extend_from_slice(&sequence[from..]);
    Some((out, moved))
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

/// An edit at its position in the edited contig, with what is needed to store
/// it as a repeat.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EditedSnp {
    pub edit: Edit,
    /// Identifies the edit across contigs: its context and change (0 without
    /// a context).
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

/// The edits of one contig (at their positions in `edited`, the edited
/// contig, see [apply_edits]) with their contexts there.
pub fn with_contexts(
    edited: &[u8],
    edits: Vec<Edit>,
) -> Vec<EditedSnp> {
    let contexts: Vec<Option<[u8; CONTEXT]>> = edits.iter().map(|e| context_at(edited, e.pos as usize)).collect();
    let present: Vec<[u8; CONTEXT]> = contexts.iter().flatten().copied().collect();
    let found = find_contexts(edited, &present);
    let mut counts = found.iter().map(|f| f.len());
    edits.into_iter().zip(contexts).map(|(edit, context)| {
        let unique = context.is_some() && counts.next() == Some(1);
        let key = context.map_or(0, |context| {
            let mut hasher = DefaultHasher::new();
            context.hash(&mut hasher);
            edit.change.hash(&mut hasher);
            hasher.finish()
        });
        EditedSnp { edit, key, context, unique }
    }).collect()
}

/// A new SNP's change as stored.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Stored {
    /// The original base's distance from the index's base, 1..=3.
    Distance(u8),
    /// A deletion's bases, put back at the position.
    Deleted(Vec<u8>),
    /// An insertion's length, removed from the position on.
    Inserted(u8),
}

impl Stored {
    fn of(change: &Change) -> Self {
        match change {
            Change::Sub { base, index_base } => Stored::Distance(base_symbol(*base, *index_base)),
            Change::Del { bases } => Stored::Deleted(bases.clone()),
            Change::Ins { bases } => Stored::Inserted(bases.len() as u8),
        }
    }
}

/// SNPs of one colour.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
struct ColourSnps {
    /// New SNPs: position and the change.
    new: Vec<(u32, Stored)>,
    /// Colours back to the reference (0: none).
    reference: u32,
    /// Indices of the shared entries in the reference's id list, increasing.
    shared: Vec<u32>,
    /// Ids of the other repeated SNPs, increasing.
    direct: Vec<u32>,
    /// Ids of all repeated SNPs: the shared entries', then the direct ones.
    repeats: Vec<u32>,
}

/// Edited SNPs of every colour.
///
/// Serialized as the arithmetic-coded stream of [SnpEdits::to_bytes].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnpEdits {
    colours: Vec<ColourSnps>,
    /// Colour and index in its new SNPs of every id.
    ids: Vec<(u32, u32)>,
}

/// Distance of `base` from `index_base`, 1..=3.
fn base_symbol(
    base: u8,
    index_base: u8,
) -> u8 {
    (base_code(base).unwrap() + 4 - base_code(index_base).unwrap()) % 4
}

fn bit_length(
    x: u64,
) -> usize {
    (64 - x.leading_zeros()) as usize
}

/// Adaptive models of the SNP stream, see the module docs.
struct Models {
    n_new: [GammaModel; 2],
    has_reference: [u32; 2],
    distance: GammaModel,
    flags: [u32; 16],
    gaps: Vec<GammaModel>,
    bases: [u32; 4],
    is_deletion: u32,
    indel_lengths: [[u32; 4]; 2],
    deleted: [u32; 4],
    n_direct: [GammaModel; 2],
    ids: Vec<GammaModel>,
}

/// Contexts of the gaps and the direct ids: the bit length of the previous
/// value (0 for none), capped.
const N_LENGTH_CONTEXTS: usize = 24;

impl Models {
    fn new() -> Self {
        Models {
            n_new: [GammaModel::new(), GammaModel::new()],
            has_reference: [C0; 2],
            distance: GammaModel::new(),
            flags: [C0; 16],
            gaps: (0..N_LENGTH_CONTEXTS).map(|_| GammaModel::new()).collect(),
            bases: [C0; 4],
            is_deletion: C0,
            indel_lengths: [[C0; 4]; 2],
            deleted: [C0; 4],
            n_direct: [GammaModel::new(), GammaModel::new()],
            ids: (0..N_LENGTH_CONTEXTS).map(|_| GammaModel::new()).collect(),
        }
    }
}

/// Codes (`colours` Some) or decodes the SNPs of every colour; returns them
/// with the id of every new SNP.
fn code_all<C: Coder>(
    c: &mut C,
    colours: Option<&[ColourSnps]>,
    overrun: &dyn Fn(&C) -> bool,
) -> Result<(Vec<ColourSnps>, Vec<(u32, u32)>), SnpDecodeError> {
    let truncated = || SnpDecodeError("SNP stream truncated".to_string());
    let mut m = Models::new();
    let n_colors = GammaModel::new().code(c, colours.map_or(0, |x| x.len()) as u64 + 1) - 1;
    let n_colors: usize = n_colors.try_into().map_err(|_| truncated())?;
    let mut out: Vec<ColourSnps> = Vec::with_capacity(n_colors.min(1 << 20));
    let mut lists: Vec<Vec<u32>> = Vec::with_capacity(n_colors.min(1 << 20));
    let mut ids: Vec<(u32, u32)> = Vec::new();
    let mut previous_empty = 1;
    for col in 0..n_colors {
        if overrun(c) {
            return Err(truncated())
        }
        let truth = colours.map(|x| &x[col]);
        let error = |what: &str| SnpDecodeError(format!("SNPs of colour {}: {}", col, what));
        let n_new = m.n_new[previous_empty].code(c, truth.map_or(0, |t| t.new.len()) as u64 + 1) - 1;
        let n_new: u32 = n_new.try_into().map_err(|_| error("too many SNPs"))?;
        let empty = (n_new == 0) as usize;
        previous_empty = empty;

        // New SNPs: the gap from the base after the previous one, and the base.
        let mut new: Vec<(u32, Stored)> = Vec::with_capacity((n_new as usize).min(1 << 20));
        let mut list: Vec<u32> = Vec::with_capacity((n_new as usize).min(1 << 20));
        let (mut next, mut previous_length) = (0_u64, 0_usize);
        for i in 0..n_new {
            let (pos, stored) = truth.map_or((0, &Stored::Distance(1)), |t| (t.new[i as usize].0, &t.new[i as usize].1));
            let gap = m.gaps[previous_length.min(N_LENGTH_CONTEXTS - 1)].code(c, pos as u64 - next + 1) - 1;
            previous_length = bit_length(gap + 1);
            let pos: u32 = (next + gap).try_into().map_err(|_| error("position beyond 32 bits"))?;
            let symbol = match stored { Stored::Distance(d) => *d as u32 - 1, _ => 3 };
            let symbol = code_tree(c, &mut m.bases, 2, symbol);
            let stored = if symbol < 3 {
                Stored::Distance(symbol as u8 + 1)
            } else {
                let is_deletion = code_bit(c, &mut m.is_deletion, matches!(stored, Stored::Deleted(_)) as u32);
                let len = match stored { Stored::Deleted(bases) => bases.len() as u32, Stored::Inserted(len) => *len as u32, _ => 1 };
                let len = code_tree(c, &mut m.indel_lengths[is_deletion as usize], 2, len - 1) + 1;
                if is_deletion == 1 {
                    let bases = (0..len as usize).map(|j| {
                        let truth = match stored { Stored::Deleted(bases) => base_code(bases[j]).unwrap() as u32, _ => 0 };
                        b"ACGT"[code_tree(c, &mut m.deleted, 2, truth) as usize]
                    }).collect();
                    Stored::Deleted(bases)
                } else {
                    Stored::Inserted(len as u8)
                }
            };
            new.push((pos, stored));
            list.push(ids.len() as u32);
            ids.push((col as u32, i));
            next = pos as u64 + 1;
            if i % 1024 == 1023 && overrun(c) {
                return Err(truncated())
            }
        }
        let ids_before = (ids.len() - n_new as usize) as u64;

        // The reference and its shared entries: a flag per entry of its id list.
        let has_reference = code_bit(c, &mut m.has_reference[empty], truth.map_or(0, |t| (t.reference > 0) as u32));
        let (mut reference, mut shared) = (0_u32, Vec::new());
        let mut repeats: Vec<u32> = Vec::new();
        if has_reference == 1 {
            let distance = m.distance.code(c, truth.map_or(1, |t| t.reference as u64));
            if distance == 0 || distance > col as u64 {
                return Err(error("reference before the first colour"))
            }
            reference = distance as u32;
            let r = col - distance as usize;
            let n_new_reference = out[r].new.len();
            let entries = &lists[r];
            let wanted = truth.map(|t| &t.shared);
            let (mut f1, mut f2, mut j) = (0, 0, 0);
            for (i, &id) in entries.iter().enumerate() {
                let flag = wanted.is_some_and(|w| j < w.len() && w[j] as usize == i);
                if flag { j += 1; }
                let context = ((f1 * 2 + f2) * 2 + (i < n_new_reference) as usize) * 2 + (i + 1 == entries.len()) as usize;
                let flag = code_bit(c, &mut m.flags[context], flag as u32) as usize;
                if flag == 1 {
                    shared.push(i as u32);
                    repeats.push(id);
                }
                (f2, f1) = (f1, flag);
            }
        }

        // Direct repeats: ids, the first as the distance back from the newest
        // earlier id, then the differences.
        let n_direct = m.n_direct[empty].code(c, truth.map_or(0, |t| t.direct.len()) as u64 + 1) - 1;
        if n_direct > ids_before {
            return Err(error("more repeats than earlier SNPs"))
        }
        let mut direct: Vec<u32> = Vec::with_capacity(n_direct as usize);
        let (mut previous, mut previous_length) = (None::<u64>, 0_usize);
        for i in 0..n_direct as usize {
            let id = truth.map_or(0, |t| t.direct[i] as u64);
            let (context, value) = match previous {
                None => (0, ids_before - id),
                Some(p) => (1 + previous_length.min(N_LENGTH_CONTEXTS - 2), id - p),
            };
            let value = m.ids[context].code(c, value);
            previous_length = bit_length(value);
            let id = match previous {
                None => ids_before.checked_sub(value),
                Some(p) => p.checked_add(value),
            }.filter(|&id| id < ids_before && value > 0).ok_or_else(|| error("invalid repeat id"))?;
            direct.push(id as u32);
            previous = Some(id);
        }
        repeats.extend_from_slice(&direct);
        list.extend_from_slice(&repeats);
        lists.push(list);
        out.push(ColourSnps { new, reference, shared, direct, repeats });
    }
    Ok((out, ids))
}

impl SnpEdits {
    /// Codes the SNPs of every colour, each sorted by position, all as new.
    pub fn new(
        snps: &[Vec<Snp>],
    ) -> Self {
        let edited: Vec<Vec<EditedSnp>> = snps.iter().map(|colour| colour.iter().map(|&snp| EditedSnp { edit: snp.into(), key: 0, context: None, unique: false }).collect()).collect();
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
        let mut colours: Vec<ColourSnps> = Vec::with_capacity(snps.len());
        let mut ids: Vec<(u32, u32)> = Vec::new();
        for (c, colour) in snps.iter().enumerate() {
            let c = c as u32;
            let candidate: Vec<Option<u32>> = colour.iter().map(|s| {
                if !s.unique { return None }
                let (id, context) = first.get(&s.key)?;
                (s.context.as_ref() == Some(context)).then_some(*id)
            }).collect();
            let mut wanted: HashSet<u32> = candidate.iter().flatten().copied().collect();
            // The reference: the earlier colour holding the most candidates.
            let mut votes: HashMap<u32, usize> = HashMap::new();
            for id in &wanted {
                for &r in holders.get(id).map_or(&[][..], |v| v.as_slice()) { *votes.entry(r).or_default() += 1; }
            }
            let reference = votes.into_iter().max_by_key(|&(r, v)| (v, r)).map(|(r, _)| r);
            let mut shared: Vec<u32> = Vec::new();
            let mut repeats: Vec<u32> = Vec::new();
            if let Some(r) = reference {
                for (i, id) in lists[r as usize].iter().enumerate() {
                    if wanted.remove(id) { shared.push(i as u32); repeats.push(*id); }
                }
            }
            // The candidates the reference does not hold are direct repeats.
            let mut direct: Vec<u32> = wanted.into_iter().collect();
            direct.sort_unstable();
            repeats.extend_from_slice(&direct);
            let repeated: HashSet<u32> = repeats.iter().copied().collect();
            let mut list: Vec<u32> = Vec::new();
            let mut new: Vec<(u32, Stored)> = Vec::new();
            for (s, cand) in colour.iter().zip(&candidate) {
                if cand.is_some_and(|id| repeated.contains(&id)) { continue }
                let id = ids.len() as u32;
                if let (Some(context), false) = (s.context, first.contains_key(&s.key)) {
                    first.insert(s.key, (id, context));
                }
                list.push(id);
                ids.push((c, new.len() as u32));
                new.push((s.edit.pos, Stored::of(&s.edit.change)));
            }
            list.extend_from_slice(&repeats);
            for &id in &list { holders.entry(id).or_default().push(c); }
            lists.push(list);
            let reference = if shared.is_empty() { 0 } else { c - reference.unwrap() };
            colours.push(ColourSnps { new, reference, shared, direct, repeats });
        }
        SnpEdits { colours, ids }
    }

    /// The arithmetic-coded stream (see the module docs).
    pub fn to_bytes(
        &self,
    ) -> Vec<u8> {
        let mut encoder = Encoder::new();
        code_all(&mut encoder, Some(&self.colours), &|_| false).expect("valid SNPs");
        encoder.finish()
    }

    /// Decodes the stream written by [SnpEdits::to_bytes].
    pub fn from_bytes(
        bytes: &[u8],
    ) -> Result<Self, SnpDecodeError> {
        let mut decoder = Decoder::new(bytes);
        let (colours, ids) = code_all(&mut decoder, None, &|d: &Decoder| d.overrun())?;
        if decoder.overrun() {
            return Err(SnpDecodeError("SNP stream truncated".to_string()))
        }
        Ok(SnpEdits { colours, ids })
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

    /// Puts the original contig of `color` back from its decoded (edited)
    /// one. `context` gives the context at a position of an earlier colour's
    /// edited contig (for the pairs in [SnpEdits::context_sources]).
    pub fn revert(
        &self,
        color: u32,
        sequence: &mut Vec<u8>,
        context: &dyn Fn(u32, u32) -> Option<[u8; CONTEXT]>,
    ) -> Result<(), SnpDecodeError> {
        let snps = self.colours.get(color as usize).ok_or_else(|| SnpDecodeError(format!("no SNP stream for colour {}", color)))?;
        // Repeats are located in the edited sequence before anything changes.
        let mut targets: Vec<[u8; CONTEXT]> = Vec::with_capacity(snps.repeats.len());
        let mut stored: Vec<&Stored> = Vec::with_capacity(snps.repeats.len());
        for &id in &snps.repeats {
            let (f, i) = self.ids[id as usize];
            let (pos, change) = &self.colours[f as usize].new[i as usize];
            targets.push(context(f, *pos).ok_or_else(|| SnpDecodeError(format!("no context at {} of colour {} for colour {}", pos, f, color)))?);
            stored.push(change);
        }
        let found = find_contexts(sequence, &targets);
        let mut changes: Vec<(usize, &Stored)> = Vec::with_capacity(snps.new.len() + targets.len());
        for (j, starts) in found.iter().enumerate() {
            if starts.len() != 1 {
                return Err(SnpDecodeError(format!("repeated SNP {} of colour {} found {} times", j, color, starts.len())))
            }
            changes.push((starts[0] + FLANK, stored[j]));
        }
        changes.extend(snps.new.iter().map(|(pos, change)| (*pos as usize, change)));
        // From the end, so that indels do not move the positions still to change.
        changes.sort_by_key(|&(pos, _)| std::cmp::Reverse(pos));
        if changes.windows(2).any(|w| w[0].0 == w[1].0) {
            return Err(SnpDecodeError(format!("two SNPs of colour {} at the same position", color)))
        }
        let invalid = |pos: usize| SnpDecodeError(format!("SNP of colour {} at {} does not fit the contig", color, pos));
        for (pos, change) in changes {
            match change {
                Stored::Distance(distance) => {
                    let index_base = sequence.get(pos).and_then(|&b| base_code(b)).ok_or_else(|| invalid(pos))?;
                    sequence[pos] = b"ACGT"[(index_base as usize + *distance as usize) % 4];
                },
                Stored::Deleted(bases) => {
                    if pos > sequence.len() {
                        return Err(invalid(pos))
                    }
                    sequence.splice(pos..pos, bases.iter().copied());
                },
                Stored::Inserted(len) => {
                    if pos + *len as usize > sequence.len() {
                        return Err(invalid(pos))
                    }
                    sequence.drain(pos..pos + *len as usize);
                },
            }
        }
        Ok(())
    }
}

impl serde::Serialize for SnpEdits {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.to_bytes().serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for SnpEdits {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let bytes = Vec::<u8>::deserialize(deserializer)?;
        SnpEdits::from_bytes(&bytes).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let colours: Vec<Vec<Snp>> = (0..50).map(|c| (0..40).map(|i| snp(i * 97 + c, b'A', b"CGT"[(i % 3) as usize])).collect()).collect();
        let bytes = SnpEdits::new(&colours).to_bytes();
        assert!(SnpEdits::from_bytes(&bytes).is_ok());
        assert!(SnpEdits::from_bytes(&bytes[..bytes.len() / 2]).is_err());
        assert!(SnpEdits::from_bytes(&[]).is_err());
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
        let contigs: Vec<Vec<u8>> = colours.iter().map(|_| original.clone()).collect();
        let edited_seqs: Vec<Vec<u8>> = colours.iter().zip(&contigs).map(|(c, s)| edit(s, c)).collect();
        let called: Vec<Vec<EditedSnp>> = colours.iter().zip(&edited_seqs).map(|(c, e)| {
            with_contexts(e, c.iter().map(|&snp| snp.into()).collect())
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

    #[test]
    fn repeats_from_colours_other_than_the_reference() {
        // Colour 2 shares one SNP with colour 0 and one with colour 1: one is
        // in its reference's list, the other a direct repeat.
        let original = random_sequence(2000, 13);
        let change = |pos: usize, by: usize| snp(pos as u32, original[pos], b"ACGT"[(base_code(original[pos]).unwrap() as usize + by) % 4]);
        let (a, b) = (change(400, 1), change(1200, 2));
        let colours: Vec<Vec<Snp>> = vec![vec![a], vec![b], vec![a, b]];
        let edited_seqs: Vec<Vec<u8>> = colours.iter().map(|c| edit(&original, c)).collect();
        let called: Vec<Vec<EditedSnp>> = colours.iter().zip(&edited_seqs).map(|(c, e)| {
            with_contexts(e, c.iter().map(|&snp| snp.into()).collect())
        }).collect();
        let edits = SnpEdits::with_repeats(&called);
        assert_eq!(edits.n_repeats(), 2);
        assert_eq!(edits.colours[2].shared.len(), 1);
        assert_eq!(edits.colours[2].direct.len(), 1);
        let edits: SnpEdits = postcard::from_bytes(&postcard::to_allocvec(&edits).unwrap()).unwrap();
        let mut sources = edits.context_sources(2);
        sources.sort();
        assert_eq!(sources, vec![(0, 400), (1, 1200)]);
        let context = |f: u32, pos: u32| context_at(&edited_seqs[f as usize], pos as usize);
        for color in 0..3 {
            let mut sequence = edited_seqs[color].clone();
            edits.revert(color as u32, &mut sequence, &context).unwrap();
            assert_eq!(sequence, original, "colour {}", color);
        }
    }

    #[test]
    fn shared_indels_are_stored_once() {
        // Two contigs with the same deletion and insertion (and the second with
        // a SNP of its own): the second stores the indels as repeats.
        let original = random_sequence(2000, 17);
        let del = Edit { pos: 500, change: Change::Del { bases: original[500..502].to_vec() } };
        let ins = Edit { pos: 1200, change: Change::Ins { bases: b"GA".to_vec() } };
        let own = Edit::from(snp(1700, original[1700], b"ACGT"[(base_code(original[1700]).unwrap() as usize + 1) % 4]));
        let colours = [vec![del.clone(), ins.clone()], vec![del, ins, own]];
        let edited: Vec<(Vec<u8>, Vec<Edit>)> = colours.iter().map(|c| apply_edits(&original, c.clone()).unwrap()).collect();
        assert_eq!(edited[0].0.len(), original.len());
        let called: Vec<Vec<EditedSnp>> = edited.iter().map(|(e, moved)| with_contexts(e, moved.clone())).collect();
        let edits: SnpEdits = postcard::from_bytes(&postcard::to_allocvec(&SnpEdits::with_repeats(&called)).unwrap()).unwrap();
        assert_eq!(edits.n_snps(), 5);
        assert_eq!(edits.n_repeats(), 2);
        let context = |f: u32, pos: u32| context_at(&edited[f as usize].0, pos as usize);
        for color in 0..2 {
            let mut sequence = edited[color].0.clone();
            edits.revert(color as u32, &mut sequence, &context).unwrap();
            assert_eq!(sequence, original, "colour {}", color);
        }
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
        // SNPs are only called after matches longer than random ones.
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
    fn snps_are_called_by_walking_the_index() {
        let reference = random_sequence(4000, 11);
        let k = 63;
        let mut opts = kbo::BuildOpts::default();
        opts.k = k;
        opts.add_revcomp = true;
        opts.build_select = true;
        let (sbwt, lcs) = kbo::build(&[reference.clone()], opts);

        // An isolated SNP, three within k of each other, one too close to the
        // end to be followed by k - 1 bases, an extra base in the contig (a
        // deletion), a missing one (an insertion), and a deleted base next to
        // a SNP (substitutions alone do not get through).
        let mut contig = reference.clone();
        let mut expected = Vec::new();
        for pos in [600, 1500, 1510, 1540, 3990] {
            let base = b"ACGT"[(base_code(reference[pos]).unwrap() as usize + 1) % 4];
            contig[pos] = base;
            if pos < 3990 {
                expected.push(Edit::from(snp(pos as u32, base, reference[pos])));
            }
        }
        let other = |pos: usize, by: usize| b"ACGT"[(base_code(reference[pos]).unwrap() as usize + by) % 4];
        let (extra, near) = (other(2499, 2), other(3300, 1));
        contig[3300] = near;
        contig.insert(3301, other(3301, 2));
        contig.remove(2900);
        contig.insert(2500, extra);
        let called: Vec<Edit> = call_snps(&contig, &sbwt, &lcs).into_iter().map(|called| called.edit).collect();
        let subs: Vec<Edit> = called.iter().filter(|e| matches!(e.change, Change::Sub { .. })).cloned().collect();
        assert_eq!(subs[..4], expected[..]);
        // The indels: the extra base and the missing one, and two edits for the
        // SNP next to an extra base (which two depends on the bases).
        assert!(called.iter().any(|e| matches!(e.change, Change::Del { .. })), "{:?}", called);
        assert!(called.iter().any(|e| matches!(e.change, Change::Ins { .. })), "{:?}", called);
        assert_eq!(called.len(), expected.len() + 4, "{:?}", called);
        // Editing gives the reference but for the SNP too close to the end.
        let (edited, moved) = apply_edits(&contig, called).unwrap();
        let mut target = reference.clone();
        target[3990] = contig[contig.len() - 10];
        assert_eq!(edited, target);
        // And reverting the stored edits gives the contig back.
        let colours = vec![with_contexts(&edited, moved)];
        let edits: SnpEdits = postcard::from_bytes(&postcard::to_allocvec(&SnpEdits::with_repeats(&colours)).unwrap()).unwrap();
        let mut decoded = edited.clone();
        edits.revert(0, &mut decoded, &|_, _| None).unwrap();
        assert_eq!(decoded, contig);
        assert!(call_snps(&reference, &sbwt, &lcs).is_empty());
    }
}
