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
//! Stored nodes and chain targets of [ChainCsr](super::ChainCsr), in node
//! order (see [stored](super::stored)).
//!
//! A stored node is its colex rank in the truncated binary code (also called
//! minimal binary) for the bound `colex_bound` = largest rank + 1: with
//! l = ⌊log2 bound⌋, the first 2^(l+1) − bound values take l bits and the
//! others l + 1 bits.
//!
//! When the graph is stored with the index ([WalkIndex]), most nodes are
//! *derived* instead. The phrase of a node is its _k_-mer left-extended by
//! the suffix length of the edge leaving it. A phrase usually ends at a
//! _k_-mer with more than one predecessor while the contig is still in the
//! index; then the next node's _k_-mer is _k_ left steps from the start of
//! the phrase in the de Bruijn graph of the index, and only the predecessor
//! taken at each step with more than one is needed. For a node X whose
//! predecessor A on its chain is a full entry (not the source, suffix length
//! implicit or at least _k_), an adaptive flag (by inner node or junction and
//! by A's suffix length implicit or explicit) says whether X is derived: then
//! each branch on the walk from A's phrase start is coded as the index of the
//! predecessor taken in the colex interval of the predecessors (adaptive, by
//! the in-degree, whether it is the first step and the choice taken the last
//! time a walk passed the same branch, if any); otherwise, and for the other
//! nodes, X's rank is stored. A walk step takes one predecessor from the
//! SBWT's inverse LF step and the others from the LCS array (they share the
//! last _k_ − 1 bases), so decoding needs the index. The encoder derives a
//! node when every _k_-mer of X followed by the first _k_ bases of A's
//! phrase is in the index and no step has more than four predecessors.
//!
//! With the index, the target of every chain is coded in node order where
//! the chain ends, instead of in the targets stream. If the chain's last node
//! A is a full entry, an adaptive flag (by A's suffix length implicit or
//! explicit) says whether the target is derived: its _k_-mer is the end of
//! the walk from A's phrase start, coded as above. A derived rank that
//! belongs to a junction reached before is that junction; any other is the
//! next new junction. Otherwise the target is coded as in "Targets" (below),
//! in the node stream, and a new target's rank is stored. The sink is never
//! derived.
//!
//! ## Targets
//!
//! Junctions are numbered in breadth-first order, so the target of a chain is
//! either the next junction not reached yet or one reached before. Without
//! the index the targets are coded first, in chain order, in their own
//! stream; with it, in the node stream (above). Each target is
//! coded with the binary arithmetic coder of [crate::codecs::arith]: an adaptive
//! flag "the next new junction"; otherwise the distance d from the newest
//! junction reached so far, as the bit length of d + 1 (a 6-bit adaptive tree)
//! and the bits below its top one, the first [TARGET_MODELLED_BITS] of them
//! adaptive given the bit length and the bits above them, the rest raw.
//!
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::codecs::arith::{code_bit, code_tree, Coder, Decoder, Encoder, C0};
use crate::codecs::bits::{BitReader, BitWriter};
use crate::graph::GraphDecodeError;
use crate::walk::{IndexWalk, WalkIndex};

/// Low bits of a target distance coded with adaptive models (the others raw).
pub const TARGET_MODELLED_BITS: u32 = 8;

/// Adaptive models of the chain targets coded as in "Stored targets", and
/// the next new junction.
pub(super) struct TargetModels {
    is_new: u32,
    length_tree: Vec<u32>,
    low: Vec<u32>,
    pub(super) next_new: u32,
}

impl TargetModels {
    pub(super) fn new() -> Self {
        TargetModels { is_new: C0, length_tree: vec![C0; 64], low: vec![C0; 34 << TARGET_MODELLED_BITS], next_new: 2 }
    }

    /// Codes `target` (or decodes one, None); returns it and whether it is
    /// the next new junction. None for an invalid target.
    fn code<C: Coder>(
        &mut self,
        c: &mut C,
        target: Option<u32>,
    ) -> Option<(u32, bool)> {
        let next_new = self.next_new;
        if let Some(target) = target {
            assert!(target <= next_new, "junctions are not numbered in breadth-first order");
        }
        if code_bit(c, &mut self.is_new, target.is_some_and(|t| t == next_new) as u32) == 1 {
            self.next_new += 1;
            return Some((next_new, true))
        }
        // d + 1 (d = distance from the newest junction) has `bits` bits.
        let d = next_new.wrapping_sub(1).wrapping_sub(target.unwrap_or(0)) as u64;
        let bits = code_tree(c, &mut self.length_tree, 6, 63 - (d + 1).leading_zeros()) + 1;
        if bits > 33 { return None }
        let mut x = 1_u64;
        for (j, b) in (0..bits - 1).rev().enumerate() {
            let truth = (((d + 1) >> b) & 1) as u32;
            let bit = if (j as u32) < TARGET_MODELLED_BITS {
                code_bit(c, &mut self.low[((bits as usize) << TARGET_MODELLED_BITS) + x as usize], truth)
            } else {
                c.code(truth, 32768)
            };
            x = (x << 1) | bit as u64;
        }
        let d = x - 1;
        if d >= next_new as u64 - 1 { return None }
        Some((next_new - 1 - d as u32, false))
    }
}

/// Codes the `targets` of chains with the arithmetic coder, or decodes `n`
/// of them (`targets` None), see the module docs. Returns the targets and the
/// number of junctions they reach (the next new junction); fails if decoded
/// targets are not in breadth-first order.
pub(super) fn code_targets<C: Coder>(
    c: &mut C,
    n: usize,
    targets: Option<&[u32]>,
) -> Result<(Vec<u32>, u32), GraphDecodeError> {
    let mut models = TargetModels::new();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let (target, _) = models.code(c, targets.map(|t| t[i])).ok_or_else(|| GraphDecodeError("stored chains: targets truncated or invalid".to_string()))?;
        out.push(target);
    }
    Ok((out, models.next_new))
}

/// Largest number of predecessors of a _k_-mer (one per base).
const MAX_IN_DEGREE: usize = 4;

/// Adaptive models of the node walks.
pub(super) struct WalkModels {
    /// Derived or not: inner node or junction x implicit or explicit suffix length.
    flags: [u32; 4],
    /// Chain target derived or not: implicit or explicit suffix length.
    ends: [u32; 2],
    /// Branch choices: the choice at the same in-branch last time (none,
    /// 0..=3) x first step or not x in-degree 2..=4 (2-bit trees).
    steps: Vec<u32>,
    /// The last choice at each in-branch, by its group's first row.
    last_at: std::collections::HashMap<u32, u8>,
}

impl WalkModels {
    pub(super) fn new() -> Self {
        WalkModels { flags: [C0; 4], ends: [C0; 2], steps: vec![C0; 5 * 2 * MAX_IN_DEGREE * 4], last_at: std::collections::HashMap::new() }
    }
}

/// Codes (or decodes, `rows` None) the branches of the walk of _k_ left steps
/// from the _k_-mer `start`; when encoding, `rows[p]` is the colex rank of
/// the _k_-mer ending at p of the 2_k_ bases of the derived node followed by
/// `start`. Returns the rank of the node reached, None for an invalid walk.
fn code_walk<C: Coder>(
    c: &mut C,
    m: &mut WalkModels,
    index: &WalkIndex,
    start: &[u8],
    rows: Option<&[usize]>,
) -> Option<u32> {
    let k = index.sbwt.k();
    let mut cur = index.sbwt.search(start)?.start;
    for step in 0..k {
        let group = index.predecessors(cur)?;
        let in_degree = group.end - group.start;
        if in_degree > MAX_IN_DEGREE {
            return None
        }
        cur = group.start;
        if in_degree > 1 {
            let truth = rows.map_or(0, |rows| (rows[2 * k - 2 - step] - group.start) as u32);
            let last = m.last_at.get(&(group.start as u32)).map_or(4, |&x| x as usize);
            let context = ((last * 2 + (step > 0) as usize) * MAX_IN_DEGREE + in_degree - 1) * 4;
            let chosen = code_tree(c, &mut m.steps[context..context + 4], 2, truth) as usize;
            if chosen >= in_degree {
                return None
            }
            m.last_at.insert(group.start as u32, chosen as u8);
            cur += chosen;
        }
    }
    u32::try_from(cur).ok()
}

/// A branch of a walk: its step, the in-degree, the first row of the
/// predecessors' group and the predecessor taken.
pub(super) type Branch = (u16, u8, u32, u8);

/// The branches of the walk of _k_ left steps through `rows` (see
/// [code_walk]; the walk starts at `rows[2k - 1]`, the phrase start), which
/// reaches the colex rank returned.
fn walk_branches(
    index: &WalkIndex,
    rows: &[usize],
) -> (Vec<Branch>, u32) {
    let k = index.sbwt.k();
    let mut cur = rows[2 * k - 1];
    let mut branches = Vec::new();
    for step in 0..k {
        let group = index.predecessors(cur).expect("derivable walk");
        let in_degree = group.end - group.start;
        cur = group.start;
        if in_degree > 1 {
            let chosen = rows[2 * k - 2 - step] - group.start;
            branches.push((step as u16, in_degree as u8, group.start as u32, chosen as u8));
            cur += chosen;
        }
    }
    (branches, cur as u32)
}

/// Codes the `branches` of a walk, as [code_walk] does.
fn code_branches<C: Coder>(
    c: &mut C,
    m: &mut WalkModels,
    branches: &[Branch],
) {
    for &(step, in_degree, group_start, chosen) in branches {
        let last = m.last_at.get(&group_start).map_or(4, |&x| x as usize);
        let context = ((last * 2 + (step > 0) as usize) * MAX_IN_DEGREE + in_degree as usize - 1) * 4;
        code_tree(c, &mut m.steps[context..context + 4], 2, chosen as u32);
        m.last_at.insert(group_start, chosen);
    }
}

/// What the encoder codes for a node or a chain's target, prepared
/// independently of the other nodes (so in parallel): None without a phrase
/// start (the rank is stored), else the predecessor's suffix length and, if
/// the node is derived, the branches of its walk.
pub(super) type Prepared = Option<(u32, Option<Vec<Branch>>)>;

/// The walk index of the encoder and its streaming index (for the matching
/// statistics that check a walk).
pub(super) type WriterWalk<'a> = (WalkIndex<'a>, StreamingIndex<'a, SbwtIndex<SubsetMatrix>, LcsArray>);

/// Encoder side of the stored nodes: the ranks, the walks and their models.
pub(super) struct NodeWriter<'a> {
    pub(super) walk: Option<WriterWalk<'a>>,
    pub(super) ranks: BitWriter,
    pub(super) walks: Encoder,
    pub(super) models: WalkModels,
    pub(super) bound: u32,
}

impl<'a> NodeWriter<'a> {
    /// The walk index and the phrase start of the predecessor `pred` (rank,
    /// suffix length; None for the source) with its suffix length, if a walk
    /// can start there. Borrows only the `walk` field, so the caller can code
    /// with the other fields meanwhile.
    fn start<'w>(
        walk: &'w Option<WriterWalk<'a>>,
        pred: Option<(u32, u32)>,
    ) -> Option<(&'w WriterWalk<'a>, Vec<u8>, u32)> {
        let (walk, (a, len)) = (walk.as_ref()?, pred?);
        walk.0.phrase_start(a, len).map(|start| (walk, start, len))
    }

    /// The rows of the k-mers of `x` followed by `start` if the walk from
    /// `start` to `x` exists: every k-mer in the index, in-degrees at most
    /// [MAX_IN_DEGREE].
    fn walk_rows(
        (index, streaming): &WriterWalk,
        start: &[u8],
        x: u32,
    ) -> Option<Vec<usize>> {
        let k = index.sbwt.k();
        let mut bases = index.sbwt.access_kmer(x as usize);
        bases.extend_from_slice(start);
        let ms = streaming.matching_statistics(&bases);
        let rows: Vec<usize> = ms.iter().map(|m| m.1.start).collect();
        let derivable = (k - 1..2 * k).all(|p| ms[p].0 == k)
            && (k..2 * k).all(|p| { let g = index.group(rows[p - 1]); g.end - g.start <= MAX_IN_DEGREE });
        derivable.then_some(rows)
    }

    /// Prepares node `x` (None: the sink as a chain's target, never derived)
    /// whose predecessor on its chain is `pred` (rank, suffix length; None
    /// for the source), see [Prepared].
    pub(super) fn prepare(
        walk: &Option<WriterWalk<'a>>,
        pred: Option<(u32, u32)>,
        x: Option<u32>,
    ) -> Prepared {
        let (walk, start, len) = Self::start(walk, pred)?;
        let branches = x.and_then(|x| Self::walk_rows(walk, &start, x).map(|rows| {
            let (branches, reached) = walk_branches(&walk.0, &rows);
            assert_eq!(reached, x, "node walk does not reach the node");
            branches
        }));
        Some((len, branches))
    }

    /// Codes node `x` of a graph stored without nodes derived by walks: its
    /// rank.
    pub(super) fn node(
        &mut self,
        x: u32,
    ) {
        self.ranks.write_truncated_binary(x, self.bound);
    }

    /// Codes the prepared node `x` ([NodeWriter::prepare]).
    pub(super) fn node_prepared(
        &mut self,
        prepared: Prepared,
        x: u32,
        is_inner: bool,
    ) {
        let Some((len, branches)) = prepared else {
            self.ranks.write_truncated_binary(x, self.bound);
            return
        };
        let flag = &mut self.models.flags[(is_inner as usize) * 2 + (len == 0) as usize];
        code_bit(&mut self.walks, flag, branches.is_some() as u32);
        match branches {
            Some(branches) => code_branches(&mut self.walks, &mut self.models, &branches),
            None => self.ranks.write_truncated_binary(x, self.bound),
        }
    }

    /// Codes the target junction `target` (colex rank `rank`, `is_new` if
    /// this chain reaches it first) of a chain, prepared
    /// ([NodeWriter::prepare]): by a walk if one exists, else by `targets`
    /// and a new target's rank.
    pub(super) fn end_prepared(
        &mut self,
        prepared: Prepared,
        target: u32,
        rank: u32,
        is_new: bool,
        targets: &mut TargetModels,
    ) {
        if let Some((len, branches)) = prepared {
            code_bit(&mut self.walks, &mut self.models.ends[(len == 0) as usize], branches.is_some() as u32);
            if let Some(branches) = branches {
                code_branches(&mut self.walks, &mut self.models, &branches);
                if is_new {
                    targets.next_new += 1;
                }
                return
            }
        }
        targets.code(&mut self.walks, Some(target)).expect("target encode");
        if is_new {
            self.ranks.write_truncated_binary(rank, self.bound);
        }
    }
}

/// Decoder side of the stored nodes, see [NodeWriter].
pub(super) struct NodeReader<'a, 'b> {
    pub(super) walk: Option<WalkIndex<'a>>,
    pub(super) ranks: BitReader<'b>,
    pub(super) walks: Decoder<'b>,
    pub(super) models: WalkModels,
    pub(super) bound: u32,
}

impl<'a> NodeReader<'a, '_> {
    /// The walk index and the phrase start of the predecessor `pred` with its
    /// suffix length, if a walk can start there; see [NodeWriter::start].
    fn start<'w>(
        walk: &'w Option<WalkIndex<'a>>,
        pred: Option<(u32, u32)>,
    ) -> Option<(&'w WalkIndex<'a>, Vec<u8>, u32)> {
        let (index, (a, len)) = (walk.as_ref()?, pred?);
        index.phrase_start(a, len).map(|start| (index, start, len))
    }

    /// Decodes a node whose predecessor is `pred`.
    pub(super) fn node(
        &mut self,
        pred: Option<(u32, u32)>,
        is_inner: bool,
    ) -> Option<u32> {
        if let Some((index, start, len)) = Self::start(&self.walk, pred) {
            let flag = &mut self.models.flags[(is_inner as usize) * 2 + (len == 0) as usize];
            if code_bit(&mut self.walks, flag, 0) == 1 {
                return code_walk(&mut self.walks, &mut self.models, index, &start, None)
            }
        }
        self.ranks.read_truncated_binary(self.bound)
    }

    /// Decodes the target of a chain whose last node is `pred`: the
    /// junction and, if this chain reaches it first, its colex rank.
    pub(super) fn end(
        &mut self,
        pred: Option<(u32, u32)>,
        targets: &mut TargetModels,
        junction_of: &std::collections::HashMap<u32, u32>,
    ) -> Option<(u32, Option<u32>)> {
        if let Some((index, start, len)) = Self::start(&self.walk, pred) {
            if code_bit(&mut self.walks, &mut self.models.ends[(len == 0) as usize], 0) == 1 {
                let rank = code_walk(&mut self.walks, &mut self.models, index, &start, None)?;
                if let Some(&junction) = junction_of.get(&rank) {
                    return Some((junction, None))
                }
                let junction = targets.next_new;
                targets.next_new += 1;
                return Some((junction, Some(rank)))
            }
        }
        let (target, is_new) = targets.code(&mut self.walks, None)?;
        if is_new {
            return Some((target, Some(self.ranks.read_truncated_binary(self.bound)?)))
        }
        Some((target, None))
    }
}
