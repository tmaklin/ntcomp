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
//! NLZ-OG graph with linear chains collapsed into single edges.
//!
//! ## Graph structure
//!
//! Represents the graph as paths through inner nodes and junctions. An *inner
//! node* is a node with exactly one incoming and one outgoing edge. Every other
//! node including the dummy source and sink are *junction nodes*. Following
//! edges from a junction through inner nodes always ends at a junction, so the
//! graph can be stored as linear chains between junctions.
//!
//! Node weights for junction nodes are colex ranks. Node weights for the
//! collapsed linear chains between junctions are concatenated in the edge
//! weight.
//!
//! ## Colours
//!
//! All edges between inner nodes always have the same colour set, as changing the
//! colour set requires introducing an incoming edge with a different colour
//! set, which converts the node into a junction. Therefore, a chain stores the
//! colex ranks of its inner nodes and the suffix lengths of its edges, and its
//! colour set only once.
//!
//! Colours only have a choice at junctions, so the choice stream has the same
//! number of entries as for the uncollapsed graph. The stored ranks refer to
//! the order of the chains in the junction's row. This means that the edges
//! between inner nodes do not need to store the colour sets if they are stored
//! in the choice stream.
//!
//! ## Stored colex ranks
//!
//! The colex ranks of the junctions and the inner nodes are coded with the
//! truncated binary code (also called minimal binary) for the bound
//! `colex_bound` = largest stored rank + 1: with l = ⌊log2 bound⌋, the first
//! 2^(l+1) − bound values take l bits and the others l + 1 bits, so a rank
//! takes log2 of the bound bits on average instead of a whole number of bits.
//! The dummy source and sink are always junctions 0 and 1 and are not stored.
//!
//! ## Stored suffix lengths
//!
//! The suffix lengths on the edges out of the dummy source node are the total
//! contig lengths and are stored as they are.
//!
//! Every other suffix length is coded in chain order with the binary
//! arithmetic coder of [crate::coder::arith], with adaptive models per
//! context:
//! - a flag: the length is implicit (0, see
//!   [SuffixLen::Implicit](super::SuffixLen));
//! - a flag: the length is above D;
//! - a length up to D: length − 1 through an adaptive binary tree;
//! - a length above D: length − D as its bit length (a 6-bit adaptive tree
//!   shared by all contexts) and the bits below its top one, the first four
//!   adaptive given the bit length and the bits above them, the rest raw,
//!
//! where D = max(_k_, [MIN_DIRECT_SUFFIX_LEN]) (see [max_direct_suffix_len]),
//! with _k_ taken from the SBWT when the graph is stored or read with the
//! index, so every short entry (< _k_) is coded directly. Without the index D
//! is [MIN_DIRECT_SUFFIX_LEN], the same as for any _k_ ≤ 255. D is stored, so
//! reading with another _k_ fails.
//!
//! Without the LCS array there is one context for all edges. With it, the
//! context of an edge is min(LCS(r), D), where r is the colex rank of the
//! edge's source node: a short entry's rank is the start of the SBWT interval
//! of its suffix, so its length is always larger than LCS(r), often by one.
//! The decoder takes LCS(r) from the LCS array or computes it from the SBWT
//! with [crate::lcs_from_sbwt].
//!
//! ## Stored out-degrees and chain lengths
//!
//! One arithmetic-coded stream ([crate::coder::arith]) holds the out-degree of
//! every junction except the sink, in junction order, and then the number of
//! inner nodes of every chain, in chain order. Each is an adaptive Elias-gamma
//! code ([GammaModel]: of the degree, and of the number of inner nodes + 1).
//! The chain lengths have one model per context: the out-degree of the
//! chain's source junction (capped at 8) and the length of the previous chain
//! of the same junction (in five classes, or none): the chains of a junction
//! are alternative paths through the same region.
//!
//! ## Stored targets
//!
//! Junctions are numbered in breadth-first order, so the target of a chain is
//! either the next junction not reached yet or one reached before. Each target
//! is coded with the binary arithmetic coder of [crate::coder::arith]: an adaptive
//! flag "the next new junction"; otherwise the distance d from the newest
//! junction reached so far, as the bit length of d + 1 (a 6-bit adaptive tree)
//! and the bits below its top one, the first [TARGET_MODELLED_BITS] of them
//! adaptive given the bit length and the bits above them, the rest raw.
//!
use super::{
    Csr,
    Choices,
    GraphDecodeError,
    SINK_COLEX,
    SINK_NODE,
    SOURCE_COLEX,
    SOURCE_NODE,
    take_edge,
    Walkable,
};
use super::bits::{BitReader, BitWriter};
use crate::coder::arith::{code_bit, code_tree, Coder, Decoder, Encoder, C0};

/// Smallest number of suffix lengths coded directly (not escaped): lengths
/// up to this are coded directly even when _k_ is smaller (left extensions
/// longer than _k_ code better directly than escaped).
pub const MIN_DIRECT_SUFFIX_LEN: u32 = 255;

/// Longest suffix length coded directly for _k_-mer size `k`
/// (None: unknown), see the module docs. Also the cap of the LCS contexts.
pub fn max_direct_suffix_len(
    k: Option<u32>,
) -> u32 {
    k.unwrap_or(0).max(MIN_DIRECT_SUFFIX_LEN)
}

/// LCS array value of a colex rank, see [super::LcsSource].
pub type LcsLookup<'a> = &'a (dyn Fn(u32) -> u32 + Sync);

/// Graph of linear chains between junctions, in compressed sparse row format.
///
/// Serialized as [StoredChains].
#[derive(Debug, PartialEq, Eq)]
pub struct ChainCsr {
    /// Colex ranks of the junction nodes.
    ///
    /// Junction with weight 0 is the dummy source node, weight 1 is the sink.
    pub junction_weights: Vec<u32>,

    /// Chains leaving junction node `j` are `row_ptrs[j]..row_ptrs[j + 1]`.
    pub row_ptrs: Vec<u32>,

    /// Index of the junction that each chain ends at.
    ///
    /// Junctions are numbered in breadth-first order.
    /// TODO describe breadth-first order in the module docs.
    pub targets: Vec<u32>,

    /// Inner nodes of linear chain `c` are
    /// `inner_weights[chain_ptrs[c]..chain_ptrs[c + 1]]`.
    pub chain_ptrs: Vec<u32>,

    /// Colex ranks of the inner nodes, concatenated over all chains.
    ///
    /// Colex ranks for chain `c` are stored at
    /// `inner_weights[chain_ptrs[c]..chain_ptrs[c + 1]]`.
    pub inner_weights: Vec<u32>,

    /// Suffix lengths of the edges in each chain, concatenated over all chains.
    ///
    /// A chain with `n` inner nodes has `n + 1` edges, so chain `c` starts at
    /// `chain_ptrs[c] + c`.
    pub suffix_lens: Vec<u32>,
}

/// Stored form of [ChainCsr], see the module docs.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct StoredChains {
    /// Number of junction nodes, including the dummy source and sink.
    pub n_junctions: u32,

    /// Every stored colex rank is below this (largest rank + 1).
    pub colex_bound: u32,

    /// Colex ranks of the junction nodes except the dummy source and sink,
    /// truncated binary for `colex_bound`.
    pub junction_weights: Vec<u8>,

    /// Colex ranks of the inner nodes, concatenated over all chains,
    /// truncated binary for `colex_bound`.
    pub inner_weights: Vec<u8>,

    /// Suffix lengths of the edges out of the dummy source node (the contig lengths).
    pub sequence_lens: Vec<u32>,

    /// Largest LCS context of the suffix length models; 0: one context for
    /// all. Otherwise the context is min(LCS(source rank), lcs_cap).
    pub lcs_cap: u16,

    /// D: suffix lengths up to this are coded directly, see the module docs.
    pub max_direct: u32,

    /// Other suffix lengths, arithmetic coded, see the module docs.
    pub suffix_lens: Vec<u8>,

    /// Out-degree of every junction except the sink, then the number of
    /// inner nodes of every chain, arithmetic coded, see the module docs.
    pub topology: Vec<u8>,

    /// Target of every chain, arithmetic coded, see the module docs.
    pub targets: Vec<u8>,
}

/// Low bits of a target distance coded with adaptive models (the others raw).
pub const TARGET_MODELLED_BITS: u32 = 8;

/// Codes the `targets` of chains with the arithmetic coder, or decodes `n`
/// of them (`targets` None), see the module docs. Returns the targets and the
/// number of junctions they reach (the next new junction); fails if decoded
/// targets are not in breadth-first order.
fn code_targets<C: Coder>(
    c: &mut C,
    n: usize,
    targets: Option<&[u32]>,
) -> Result<(Vec<u32>, u32), GraphDecodeError> {
    let invalid = || GraphDecodeError("stored chains: targets truncated or invalid".to_string());
    let mut is_new = C0;
    let mut length_tree = vec![C0; 64];
    let mut low = vec![C0; 34 << TARGET_MODELLED_BITS];
    let mut out = Vec::with_capacity(n);
    let mut next_new = 2_u32;
    for i in 0..n {
        let target = targets.map_or(0, |t| t[i]);
        if targets.is_some() {
            assert!(target <= next_new, "junctions are not numbered in breadth-first order");
        }
        if code_bit(c, &mut is_new, (target == next_new) as u32) == 1 {
            out.push(next_new);
            next_new += 1;
            continue
        }
        // d + 1 (d = distance from the newest junction) has `bits` bits.
        let d = next_new.wrapping_sub(1).wrapping_sub(target) as u64;
        let bits = code_tree(c, &mut length_tree, 6, 63 - (d + 1).leading_zeros()) + 1;
        if bits > 33 { return Err(invalid()) }
        let mut x = 1_u64;
        for (j, b) in (0..bits - 1).rev().enumerate() {
            let truth = (((d + 1) >> b) & 1) as u32;
            let bit = if (j as u32) < TARGET_MODELLED_BITS {
                code_bit(c, &mut low[((bits as usize) << TARGET_MODELLED_BITS) + x as usize], truth)
            } else {
                c.code(truth, 32768)
            };
            x = (x << 1) | bit as u64;
        }
        let d = x - 1;
        if d >= next_new as u64 - 1 { return Err(invalid()) }
        out.push(next_new - 1 - d as u32);
    }
    Ok((out, next_new))
}

/// Adaptive Elias-gamma code of a number w ≥ 1: the bit length of w through a
/// 6-bit adaptive tree and the bits below its top one, the first
/// [GammaModel::MODELLED_BITS] adaptive given the bit length and the bits
/// above them, the rest raw.
struct GammaModel {
    bits: Vec<u32>,
    low: Vec<u32>,
}

impl GammaModel {
    const MODELLED_BITS: u32 = 4;

    fn new() -> Self {
        GammaModel { bits: vec![C0; 64], low: vec![C0; 65 << Self::MODELLED_BITS] }
    }

    /// Codes (or decodes, `w` ignored) a number w ≥ 1.
    fn code<C: Coder>(
        &mut self,
        c: &mut C,
        w: u64,
    ) -> u64 {
        let n = code_tree(c, &mut self.bits, 6, 63 - w.max(1).leading_zeros()) + 1;
        let mut x = 1_u64;
        for (j, b) in (0..n - 1).rev().enumerate() {
            let truth = ((w >> b) & 1) as u32;
            let bit = if (j as u32) < Self::MODELLED_BITS {
                code_bit(c, &mut self.low[((n as usize) << Self::MODELLED_BITS) + x as usize], truth)
            } else {
                c.code(truth, 32768)
            };
            x = (x << 1) | bit as u64;
        }
        x
    }
}

/// Context of a chain length: the source junction's out-degree (capped at 8)
/// and the number of inner nodes of the previous chain of the junction.
fn chain_len_context(
    degree: u64,
    previous: Option<u64>,
) -> usize {
    let previous = match previous { None => 5, Some(0) => 0, Some(1) => 1, Some(2..=3) => 2, Some(4..=15) => 3, Some(_) => 4 };
    degree.min(8) as usize * 6 + previous
}

/// Codes the out-degrees of the junctions (`degrees`, 0 for the sink) and the
/// number of inner nodes of the chains (`inner`), or decodes them for
/// `n_junctions` junctions (`shape` None), see the module docs. None if decoded
/// values are out of range.
fn code_topology<C: Coder>(
    c: &mut C,
    n_junctions: usize,
    shape: Option<(&[u64], &[u64])>,
) -> Option<(Vec<u64>, Vec<u64>)> {
    let mut degree_model = GammaModel::new();
    let mut degrees = vec![0_u64; n_junctions];
    let mut n_chains = 0_u64;
    for junction in (0..n_junctions).filter(|&j| j != SINK_NODE as usize) {
        degrees[junction] = degree_model.code(c, shape.map_or(1, |s| s.0[junction]));
        n_chains += degrees[junction];
        if n_chains > u32::MAX as u64 { return None }
    }
    let mut models: Vec<GammaModel> = (0..9 * 6).map(|_| GammaModel::new()).collect();
    let mut inner = Vec::with_capacity(n_chains as usize);
    for &degree in &degrees {
        let mut previous = None;
        for _ in 0..degree {
            let w = models[chain_len_context(degree, previous)].code(c, shape.map_or(1, |s| s.1[inner.len()] + 1));
            inner.push(w - 1);
            previous = Some(w - 1);
        }
    }
    Some((degrees, inner))
}

/// Adaptive models of the suffix lengths (see the module docs) for
/// `n_contexts` contexts and direct lengths up to `max_direct`.
struct SuffixModel {
    max_direct: u32,
    /// Bits of a direct length − 1.
    bits: u32,
    implicit: Vec<u32>,
    escaped: Vec<u32>,
    trees: Vec<u32>,
    long: GammaModel,
}

impl SuffixModel {
    fn new(
        n_contexts: usize,
        max_direct: u32,
    ) -> Self {
        let bits = 32 - (max_direct - 1).leading_zeros();
        SuffixModel {
            max_direct, bits,
            implicit: vec![C0; n_contexts],
            escaped: vec![C0; n_contexts],
            trees: vec![C0; n_contexts << bits],
            long: GammaModel::new(),
        }
    }

    /// Codes (or decodes, `len` ignored) the suffix length of an edge in `context`.
    fn code<C: Coder>(
        &mut self,
        c: &mut C,
        context: usize,
        len: u32,
    ) -> Option<u32> {
        if code_bit(c, &mut self.implicit[context], (len == 0) as u32) == 1 {
            return Some(0)
        }
        if code_bit(c, &mut self.escaped[context], (len > self.max_direct) as u32) == 0 {
            let tree = &mut self.trees[context << self.bits..(context + 1) << self.bits];
            let value = code_tree(c, tree, self.bits, len.wrapping_sub(1)) + 1;
            return (value <= self.max_direct).then_some(value)
        }
        let x = self.long.code(c, len.wrapping_sub(self.max_direct) as u64);
        u32::try_from(self.max_direct as u64).ok()?.checked_add(u32::try_from(x).ok()?)
    }
}

/// Codes `ranks`, all below `bound`, with the truncated binary code.
fn pack_colex(
    ranks: &[u32],
    bound: u32,
) -> Vec<u8> {
    let mut writer = BitWriter::new();
    for &rank in ranks {
        writer.write_truncated_binary(rank, bound);
    }
    writer.finish()
}

/// Reads `n` ranks coded by [pack_colex]; None unless `bytes` holds exactly
/// that many (up to the padding of the last byte).
fn unpack_colex(
    bytes: &[u8],
    n: usize,
    bound: u32,
) -> Option<Vec<u32>> {
    if bound == 0 {
        return None
    }
    let mut reader = BitReader::new(bytes);
    let ranks = (0..n).map(|_| reader.read_truncated_binary(bound)).collect::<Option<Vec<u32>>>()?;
    (reader.bits_read().div_ceil(8) == bytes.len()).then_some(ranks)
}

impl ChainCsr {
    /// Collapses the linear chains in `csr`.
    pub fn from_csr(
        csr: &Csr<u32>,
    ) -> Self {
        let n_nodes = csr.node_weights.len();
        let n_edges = csr.column_indices.len();

        let mut in_degree = vec![0_u32; n_nodes];
        for &target in &csr.column_indices {
            in_degree[target as usize] += 1;
        }

        let is_junction: Vec<bool> = (0..n_nodes).map(|node| {
            let (start, end) = csr.row(node);
            let is_dummy = node == SOURCE_NODE as usize || node == SINK_NODE as usize;
            is_dummy || in_degree[node] != 1 || end - start != 1
        }).collect();

        // Node index -> junction index in breadth-first order, u32::MAX until
        // reached. `order` is the queue of the search and lists the nodes of
        // the junctions in index order.
        let mut junction_index = vec![u32::MAX; n_nodes];
        junction_index[SOURCE_NODE as usize] = 0;
        junction_index[SINK_NODE as usize] = 1;
        let mut order: Vec<usize> = vec![SOURCE_NODE as usize, SINK_NODE as usize];

        let mut row_ptrs = Vec::with_capacity(n_nodes + 1);
        let mut targets = Vec::new();
        let mut chain_ptrs = vec![0_u32];
        let mut inner_weights = Vec::new();
        let mut suffix_lens = Vec::with_capacity(n_edges);

        row_ptrs.push(0);
        let mut head = 0;
        while head < order.len() {
            let node = order[head];
            head += 1;

            let (start, end) = csr.row(node);
            for pos in start..end {
                suffix_lens.push(csr.suffix_lens[pos]);

                let mut next = csr.column_indices[pos] as usize;
                while !is_junction[next] {
                    let edge = csr.row(next).0;
                    inner_weights.push(csr.node_weights[next]);
                    suffix_lens.push(csr.suffix_lens[edge]);
                    next = csr.column_indices[edge] as usize;
                }

                if junction_index[next] == u32::MAX {
                    junction_index[next] = order.len() as u32;
                    order.push(next);
                }
                targets.push(junction_index[next]);
                chain_ptrs.push(inner_weights.len() as u32);
            }
            row_ptrs.push(targets.len() as u32);
        }
        assert_eq!(order.len(), is_junction.iter().filter(|x| **x).count(), "junction not reachable from the source");

        ChainCsr {
            junction_weights: order.iter().map(|&node| csr.node_weights[node]).collect(),
            row_ptrs,
            targets,
            chain_ptrs,
            inner_weights,
            suffix_lens,
        }
    }

    /// Converts into the stored form.
    ///
    /// Panics if the junctions are not numbered in breadth-first order.
    pub fn to_stored(
        &self,
    ) -> StoredChains {
        self.to_stored_with(None, 0, None)
    }

    /// Colex rank of the source of every edge, in chain order (the junction
    /// for the first edge of a chain, else the previous inner node).
    fn edge_sources(
        &self,
    ) -> Vec<u32> {
        let mut sources = Vec::with_capacity(self.suffix_lens.len());
        for junction in 0..self.junction_weights.len() {
            let (start, end) = self.row(junction);
            for chain in start..end {
                sources.push(self.junction_weights[junction]);
                sources.extend_from_slice(self.inner_nodes(chain));
            }
        }
        sources
    }

    /// Same as [ChainCsr::to_stored], with one suffix length code per LCS
    /// value of the edge's source rank (capped at `lcs_cap`) if `lcs` gives
    /// the LCS array value of a colex rank, and direct suffix length symbols
    /// for _k_-mer size `k`, see the module docs.
    pub fn to_stored_with(
        &self,
        lcs: Option<LcsLookup>,
        lcs_cap: u32,
        k: Option<u32>,
    ) -> StoredChains {
        let lcs_cap = if lcs.is_some() { lcs_cap.max(1) } else { 0 };
        let max_direct = max_direct_suffix_len(k);
        let n_junctions = self.junction_weights.len();

        let (sink_start, sink_end) = self.row(SINK_NODE as usize);
        assert_eq!(sink_start, sink_end, "sink has outgoing chains");
        let degrees: Vec<u64> = (0..n_junctions).map(|j| { let (start, end) = self.row(j); (end - start) as u64 }).collect();
        assert!(degrees.iter().enumerate().all(|(j, &d)| d > 0 || j == SINK_NODE as usize), "a junction has no outgoing chains");
        let inner: Vec<u64> = (0..self.n_chains()).map(|chain| self.inner_nodes(chain).len() as u64).collect();
        let mut topology = Encoder::new();
        code_topology(&mut topology, n_junctions, Some((&degrees, &inner))).expect("topology encode");
        let mut targets = Encoder::new();
        let (_, next_new) = code_targets(&mut targets, self.n_chains(), Some(&self.targets)).expect("targets encode");
        assert_eq!(next_new as usize, n_junctions, "junctions are not numbered in breadth-first order");

        // Chains out of the source come first; their first edge holds the contig length.
        let n_source_chains = self.row(SOURCE_NODE as usize).1;
        let sources = if lcs_cap > 0 { self.edge_sources() } else { Vec::new() };
        let mut sequence_lens = Vec::with_capacity(n_source_chains);
        let mut model = SuffixModel::new(lcs_cap as usize + 1, max_direct);
        let mut suffix_lens = Encoder::new();
        let mut edge = 0;
        for chain in 0..self.n_chains() {
            for (i, &len) in self.chain_suffix_lens(chain).iter().enumerate() {
                let source = if lcs_cap > 0 { sources[edge] } else { 0 };
                edge += 1;
                if chain < n_source_chains && i == 0 {
                    sequence_lens.push(len);
                    continue;
                }
                let context = if lcs_cap > 0 { lcs.unwrap()(source).min(lcs_cap) as usize } else { 0 };
                model.code(&mut suffix_lens, context, len).expect("suffix length encode");
            }
        }

        // Real colex ranks only; the dummy source and sink are junctions 0 and 1.
        assert_eq!(self.junction_weights[..2], [SOURCE_COLEX, SINK_COLEX], "junctions 0 and 1 are not the source and sink");
        let real_junctions = &self.junction_weights[2..];
        let max_colex = real_junctions.iter().chain(&self.inner_weights).copied().max().unwrap_or(0);
        assert!(max_colex < SOURCE_COLEX, "colex rank collides with dummy node ranks");
        let colex_bound = max_colex + 1;

        StoredChains {
            n_junctions: n_junctions.try_into().expect("too many junctions"),
            colex_bound,
            junction_weights: pack_colex(real_junctions, colex_bound),
            inner_weights: pack_colex(&self.inner_weights, colex_bound),
            sequence_lens,
            lcs_cap: lcs_cap as u16,
            max_direct,
            suffix_lens: suffix_lens.finish(),
            topology: topology.finish(),
            targets: targets.finish(),
        }
    }

    /// Converts from the stored form.
    pub fn from_stored(
        stored: StoredChains,
    ) -> Result<Self, GraphDecodeError> {
        Self::from_stored_with(stored, None, None)
    }

    /// Same as [ChainCsr::from_stored]; `lcs` must give the LCS array value
    /// of a colex rank if the chains were stored with LCS contexts, and `k`
    /// must be the _k_-mer size they were stored with (None for _k_ ≤ 255).
    pub fn from_stored_with(
        stored: StoredChains,
        lcs: Option<LcsLookup>,
        k: Option<u32>,
    ) -> Result<Self, GraphDecodeError> {
        let max_direct = max_direct_suffix_len(k);
        if stored.max_direct != max_direct {
            return Err(GraphDecodeError(format!("stored chains: suffix lengths coded for D = {}, expected {} (stored with another k?)", stored.max_direct, max_direct)))
        }
        let truncated = |what: &str| GraphDecodeError(format!("stored chains: {} truncated or invalid", what));
        let n_junctions = stored.n_junctions as usize;
        if n_junctions < 2 {
            return Err(GraphDecodeError("stored chains: missing source or sink".to_string()))
        }
        let colex_bound = stored.colex_bound;
        if colex_bound == 0 || colex_bound > SOURCE_COLEX {
            return Err(truncated("colex rank bound"))
        }
        let mut junction_weights = vec![SOURCE_COLEX, SINK_COLEX];
        junction_weights.extend(unpack_colex(&stored.junction_weights, n_junctions - 2, colex_bound).ok_or_else(|| truncated("junction colex ranks"))?);

        let mut coded = Decoder::new(&stored.topology);
        let (degrees, inner) = code_topology(&mut coded, n_junctions, None).ok_or_else(|| truncated("out-degrees and chain lengths"))?;
        if coded.overrun() {
            return Err(truncated("out-degrees and chain lengths"))
        }
        let mut row_ptrs = Vec::with_capacity(n_junctions + 1);
        row_ptrs.push(0_u32);
        for degree in degrees {
            let ptr = (*row_ptrs.last().unwrap() as u64 + degree).try_into().map_err(|_| truncated("out-degrees"))?;
            row_ptrs.push(ptr);
        }
        let n_chains = *row_ptrs.last().unwrap() as usize;

        let mut chain_ptrs = Vec::with_capacity(n_chains + 1);
        chain_ptrs.push(0_u32);
        for len in inner {
            let ptr = (*chain_ptrs.last().unwrap() as u64 + len).try_into().map_err(|_| truncated("chain lengths"))?;
            chain_ptrs.push(ptr);
        }
        let n_inner = *chain_ptrs.last().unwrap() as usize;
        let inner_weights = unpack_colex(&stored.inner_weights, n_inner, colex_bound).ok_or_else(|| truncated("inner colex ranks"))?;

        let n_source_chains = row_ptrs[1] as usize;
        if stored.sequence_lens.len() != n_source_chains {
            return Err(truncated("contig lengths"))
        }
        let lcs_cap = stored.lcs_cap as u32;
        let lcs = match (lcs_cap, lcs) {
            (0, _) => None,
            (_, Some(lcs)) => Some(lcs),
            (_, None) => return Err(GraphDecodeError("stored chains: the suffix lengths need the LCS array".to_string())),
        };
        let mut model = SuffixModel::new(lcs_cap as usize + 1, max_direct);
        let mut coded = Decoder::new(&stored.suffix_lens);
        let mut suffix_lens = Vec::with_capacity(n_inner + n_chains);
        let mut junction = 0;
        for chain in 0..n_chains {
            while row_ptrs[junction + 1] as usize <= chain {
                junction += 1;
            }
            let n_edges = (chain_ptrs[chain + 1] - chain_ptrs[chain]) as usize + 1;
            for i in 0..n_edges {
                if chain < n_source_chains && i == 0 {
                    suffix_lens.push(stored.sequence_lens[chain]);
                    continue;
                }
                let context = match lcs {
                    None => 0,
                    Some(lcs) => {
                        let source = if i == 0 { junction_weights[junction] } else { inner_weights[chain_ptrs[chain] as usize + i - 1] };
                        lcs(source).min(lcs_cap) as usize
                    },
                };
                suffix_lens.push(model.code(&mut coded, context, 0).ok_or_else(|| truncated("suffix lengths"))?);
            }
        }

        if coded.overrun() {
            return Err(truncated("suffix lengths"))
        }

        let mut coded = Decoder::new(&stored.targets);
        let (targets, next_new) = code_targets(&mut coded, n_chains, None)?;
        if next_new as usize != n_junctions || coded.overrun() {
            return Err(truncated("targets"))
        }

        Ok(ChainCsr {
            junction_weights,
            row_ptrs,
            targets,
            chain_ptrs,
            inner_weights,
            suffix_lens,
        })
    }

    /// Number of linear chains.
    pub fn n_chains(&self) -> usize {
        self.targets.len()
    }

    /// Range of chains leaving junction `junction`.
    pub fn row(
        &self,
        junction: usize,
    ) -> (usize, usize) {
        (self.row_ptrs[junction] as usize, self.row_ptrs[junction + 1] as usize)
    }

    /// Colex ranks of the inner nodes of `chain`.
    pub fn inner_nodes(
        &self,
        chain: usize,
    ) -> &[u32] {
        &self.inner_weights[self.chain_ptrs[chain] as usize..self.chain_ptrs[chain + 1] as usize]
    }

    /// Suffix lengths of the edges in `chain`, one more than its inner nodes.
    pub fn chain_suffix_lens(
        &self,
        chain: usize,
    ) -> &[u32] {
        let start = self.chain_ptrs[chain] as usize + chain;
        let end = self.chain_ptrs[chain + 1] as usize + chain + 1;
        &self.suffix_lens[start..end]
    }

    /// Expands the chains back into an uncollapsed graph.
    ///
    /// Junctions keep their indices and inner nodes are numbered after them.
    /// Edges out of a junction keep the order of the chains, so the choice
    /// stream stays valid for the expanded graph.
    pub fn to_csr(
        &self,
    ) -> Csr<u32> {
        let n_junctions = self.junction_weights.len() as u32;
        let mut node_weights = self.junction_weights.clone();
        node_weights.extend(self.inner_weights.iter().copied());

        // (source, target, suffix length); junction rows come first in chain order.
        let mut edges: Vec<(u32, u32, u32)> = Vec::with_capacity(self.suffix_lens.len());
        for junction in 0..n_junctions as usize {
            let (start, end) = self.row(junction);
            for chain in start..end {
                let suffix_lens = self.chain_suffix_lens(chain);
                let mut prev = junction as u32;
                for (i, suffix_len) in suffix_lens.iter().enumerate().take(self.inner_nodes(chain).len()) {
                    let node = n_junctions + self.chain_ptrs[chain] + i as u32;
                    edges.push((prev, node, *suffix_len));
                    prev = node;
                }
                edges.push((prev, self.targets[chain], *suffix_lens.last().unwrap()));
            }
        }
        edges.sort_by_key(|edge| edge.0);

        let mut row_ptrs = vec![0_u32; node_weights.len() + 1];
        for edge in &edges {
            row_ptrs[edge.0 as usize + 1] += 1;
        }
        for i in 0..node_weights.len() {
            row_ptrs[i + 1] += row_ptrs[i];
        }

        Csr {
            node_weights,
            row_ptrs,
            column_indices: edges.iter().map(|edge| edge.1).collect(),
            suffix_lens: edges.iter().map(|edge| edge.2).collect(),
        }
    }
}

impl Walkable for ChainCsr {
    /// Walks the path of `color` from the source to the sink.
    ///
    /// Same as the walk of the uncollapsed graph; the bound from [Walkable]
    /// counts uncollapsed edges.
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError> {
        let mut reader = choices
            .reader(color)
            .ok_or_else(|| GraphDecodeError(format!("no choice stream for colour {}", color)))?;
        let mut junction = SOURCE_NODE as usize;
        let mut sequence_length = 0_u32;
        let mut dictionary: Vec<(u32, u32)> = Vec::new();
        let mut n_edges = 0_usize;

        loop {
            let chain = take_edge(self.row(junction), &mut reader, color, junction)?;
            let inner = self.inner_nodes(chain);
            let suffix_lens = self.chain_suffix_lens(chain);

            if junction == SOURCE_NODE as usize {
                sequence_length = suffix_lens[0];
            }
            n_edges += suffix_lens.len();
            let max_edges = (sequence_length as usize).saturating_add(1);
            if n_edges > max_edges {
                return Err(GraphDecodeError(format!("colour {} did not reach the sink in {} edges", color, max_edges)))
            }

            if junction != SOURCE_NODE as usize {
                dictionary.push((self.junction_weights[junction], suffix_lens[0]));
            }
            dictionary.extend(inner.iter().copied().zip(suffix_lens[1..].iter().copied()));

            junction = self.targets[chain] as usize;
            if junction == SINK_NODE as usize {
                if !reader.is_finished() {
                    return Err(GraphDecodeError(format!("unused choices left for colour {}", color)))
                }
                return Ok((sequence_length, dictionary))
            }
        }
    }
}

impl serde::Serialize for ChainCsr {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.to_stored().serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for ChainCsr {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let stored = StoredChains::deserialize(deserializer)?;
        ChainCsr::from_stored(stored).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::super::*;
    use super::{max_direct_suffix_len, MIN_DIRECT_SUFFIX_LEN};

    use core::ops::Range;

    fn expected(dict: &[(usize, Range<usize>)]) -> Vec<(u32, u32)> {
        dict.iter().map(|(s, r)| (r.start as u32, *s as u32)).collect()
    }

    // Shared linear stretches (50 -> 60 -> 70) and a repeat at node 10.
    fn dicts() -> Vec<Vec<(usize, Range<usize>)>> {
        vec![
            vec![(5, 10..11), (1, 50..51), (2, 60..61), (3, 70..71), (5, 10..11), (4, 30..31), (5, 10..11), (3, 20..21)],
            vec![(5, 10..11), (4, 30..31), (5, 10..11), (1, 50..51), (2, 60..61), (3, 70..71), (2, 40..41)],
            vec![(6, 80..81), (1, 50..51), (2, 60..61), (3, 70..71), (2, 40..41), (7, 90..91)],
            vec![(9, 95..96)],
        ]
    }

    fn builder(dicts: &[Vec<(usize, Range<usize>)>]) -> GraphBuilder {
        let mut builder = GraphBuilder::new();
        for d in dicts {
            builder.add_path(d.iter().map(|x| x.0).sum(), d);
        }
        builder
    }

    #[test]
    fn chain_walk_recovers_paths() {
        let dicts = dicts();
        let (chains, choices) = builder(&dicts).to_chains();
        let (csr, plain_choices) = builder(&dicts).to_csr();

        assert!(chains.n_chains() < csr.column_indices.len());
        assert_eq!(chains.suffix_lens.len(), csr.suffix_lens.len());
        assert_eq!(choices.n_choices(), plain_choices.n_choices());

        let bytes = postcard::to_allocvec(&chains).unwrap();
        let chains: ChainCsr = postcard::from_bytes(&bytes).unwrap();

        for (color, d) in dicts.iter().enumerate().rev() {
            let (len, decoded) = chains.walk_path(&choices, color as u32).unwrap();
            assert_eq!(decoded, expected(d));
            assert_eq!(len as usize, d.iter().map(|x| x.0).sum::<usize>());
        }
    }

    #[test]
    fn expanded_chains_decode_with_the_same_choices() {
        let dicts = dicts();
        let (chains, choices) = builder(&dicts).to_chains();
        let expanded = chains.to_csr();

        for (color, d) in dicts.iter().enumerate() {
            let (_, decoded) = expanded.walk_path(&choices, color as u32).unwrap();
            assert_eq!(decoded, expected(d));
        }
    }

    #[test]
    fn both_layouts_have_the_same_choices() {
        let dicts = dicts();
        let (_, chain_choices) = builder(&dicts).to_chains();
        let (_, plain_choices) = builder(&dicts).to_csr();
        assert_eq!(chain_choices.streams, plain_choices.streams);
    }

    #[test]
    fn colorsets_are_recovered_from_the_walks() {
        let dicts = dicts();
        // (source colex, target colex, suffix length) -> colours, with the dummy
        // nodes' colex ranks at the ends of every path.
        let mut expected_sets: std::collections::HashMap<(u32, u32, u32), Vec<u32>> = std::collections::HashMap::new();
        for (color, d) in dicts.iter().enumerate() {
            let mut nodes = vec![SOURCE_COLEX];
            nodes.extend(d.iter().map(|x| x.1.start as u32));
            nodes.push(SINK_COLEX);
            let mut weights = vec![d.iter().map(|x| x.0 as u32).sum::<u32>()];
            weights.extend(d.iter().map(|x| x.0 as u32));
            for (i, w) in weights.iter().enumerate() {
                let colors = expected_sets.entry((nodes[i], nodes[i + 1], *w)).or_default();
                if colors.last() != Some(&(color as u32)) {
                    colors.push(color as u32);
                }
            }
        }

        let (csr, choices) = builder(&dicts).to_csr();
        let (chains, chain_choices) = builder(&dicts).to_chains();
        for (csr, choices) in [(csr, choices), (chains.to_csr(), chain_choices)] {
            let colorsets = csr.colorsets(&choices, dicts.len() as u32).unwrap();
            assert_eq!(colorsets.len(), expected_sets.len());
            for node in 0..csr.node_weights.len() {
                let (start, end) = csr.row(node);
                for pos in start..end {
                    let key = (csr.node_weights[node], csr.node_weights[csr.column_indices[pos] as usize], csr.suffix_lens[pos]);
                    assert_eq!(colorsets[pos], expected_sets[&key]);
                }
            }
        }
    }

    #[test]
    fn stored_chains_roundtrip() {
        let dicts = dicts();
        let (chains, _) = builder(&dicts).to_chains();
        assert_eq!(chains.junction_weights[..2], [SOURCE_COLEX, SINK_COLEX]);

        let stored = chains.to_stored();
        assert_eq!(ChainCsr::from_stored(stored).unwrap(), chains);

        for field in 0..9 {
            let mut stored = chains.to_stored();
            match field {
                0 => stored.topology.clear(),
                1 => { stored.topology.pop(); },
                2 => stored.suffix_lens.clear(),
                3 => { stored.sequence_lens.pop(); },
                4 => { stored.junction_weights.pop(); },
                5 => { stored.inner_weights.pop(); },
                6 => stored.colex_bound = 0,
                7 => stored.colex_bound = 1,
                _ => stored.targets.clear(),
            }
            assert!(ChainCsr::from_stored(stored).is_err());
        }
    }

    #[test]
    fn colex_ranks_use_the_bits_they_need() {
        // Largest rank 95: bound 96, so ranks below 32 take 6 bits and the others 7.
        let (chains, _) = builder(&dicts()).to_chains();
        let stored = chains.to_stored();
        assert_eq!(stored.colex_bound, 96);
        assert_eq!(stored.n_junctions as usize, chains.junction_weights.len());
        let bits = |ranks: &[u32]| ranks.iter().map(|&r| if r < 32 { 6 } else { 7 }).sum::<usize>().div_ceil(8);
        assert_eq!(stored.junction_weights.len(), bits(&chains.junction_weights[2..]));
        assert_eq!(stored.inner_weights.len(), bits(&chains.inner_weights));

        // Ranks of 0 only (no bits), a power-of-two bound, and ranks just below the dummy ranks.
        for (ranks, bound) in [(vec![0, 0], 1), (vec![5, 63, 0], 64), (vec![1, SOURCE_COLEX - 1, 7], SOURCE_COLEX)] {
            let d: Vec<Vec<(usize, Range<usize>)>> = vec![ranks.iter().map(|&r| (3, r as usize..r as usize + 1)).collect()];
            let (chains, choices) = builder(&d).to_chains();
            let stored = chains.to_stored();
            assert_eq!(stored.colex_bound, bound);
            if bound == 1 {
                assert!(stored.junction_weights.is_empty() && stored.inner_weights.is_empty());
            }
            let decoded = ChainCsr::from_stored(stored).unwrap();
            assert_eq!(decoded, chains);
            let (_, path) = decoded.walk_path(&choices, 0).unwrap();
            assert_eq!(path.iter().map(|x| x.0).collect::<Vec<u32>>(), ranks);
        }
    }

    #[test]
    fn implicit_and_long_suffix_lengths_roundtrip() {
        // Lengths around the direct/long symbol boundary, some flagged implicit.
        // The edge 40 -> 50 with length 12 is used by both colours but flagged
        // only by colour 1, so it is implicit for both.
        let dicts: Vec<Vec<(usize, Range<usize>)>> = vec![
            vec![(255, 10..11), (256, 20..21), (70000, 30..31), (12, 40..41), (300, 50..51)],
            vec![(256, 20..21), (12, 40..41), (255, 10..11), (1, 60..61), (12, 40..41), (300, 50..51)],
        ];
        let flags = vec![
            vec![false, true, true, false, false],
            vec![true, false, false, false, true, false],
        ];
        let expected_implicit = vec![
            vec![false, true, true, true, false],
            vec![true, false, false, false, true, false],
        ];
        let mut builder = GraphBuilder::new();
        for (d, f) in dicts.iter().zip(&flags) {
            builder.add_path_with_flags(d.iter().map(|x| x.0).sum(), d, f);
        }
        let (chains, choices) = builder.to_chains();

        let bytes = postcard::to_allocvec(&chains).unwrap();
        let decoded: ChainCsr = postcard::from_bytes(&bytes).unwrap();
        assert_eq!(decoded, chains);

        for (color, (d, implicit)) in dicts.iter().zip(&expected_implicit).enumerate() {
            let (len, path) = decoded.walk_path(&choices, color as u32).unwrap();
            assert_eq!(len as usize, d.iter().map(|x| x.0).sum::<usize>());
            let expected: Vec<(u32, u32)> = d.iter().zip(implicit).map(|((s, r), &imp)| {
                (r.start as u32, if imp { IMPLICIT_SUFFIX_LEN } else { *s as u32 })
            }).collect();
            assert_eq!(path, expected);
        }
    }

    #[test]
    fn direct_suffix_symbols_follow_k() {
        // Lengths around 255 and around a k of 300.
        let dicts: Vec<Vec<(usize, Range<usize>)>> = vec![
            vec![(254, 10..11), (256, 20..21), (299, 30..31), (300, 40..41), (301, 50..51), (5000, 60..61)],
            vec![(256, 20..21), (300, 40..41), (12, 70..71)],
        ];
        let mut builder = GraphBuilder::new();
        for d in &dicts {
            builder.add_path(d.iter().map(|x| x.0).sum(), d);
        }
        let (chains, _) = builder.to_chains();

        // For k <= 255 the coding is the same as without k.
        for k in [Some(31), Some(255)] {
            assert_eq!(postcard::to_allocvec(&chains.to_stored_with(None, 0, k)).unwrap(), postcard::to_allocvec(&chains.to_stored()).unwrap());
        }
        assert_eq!(max_direct_suffix_len(Some(63)), MIN_DIRECT_SUFFIX_LEN);
        assert_eq!(max_direct_suffix_len(Some(300)), 300);

        // With k = 300 every length up to 300 is coded directly, and the chains
        // only decode with the same k.
        let stored = chains.to_stored_with(None, 0, Some(300));
        assert_eq!(stored.max_direct, 300);
        assert_eq!(chains.to_stored().max_direct, MIN_DIRECT_SUFFIX_LEN);
        assert_eq!(ChainCsr::from_stored_with(chains.to_stored_with(None, 0, Some(300)), None, Some(300)).unwrap(), chains);
        assert!(ChainCsr::from_stored_with(chains.to_stored_with(None, 0, Some(300)), None, None).is_err());
        assert!(ChainCsr::from_stored_with(chains.to_stored(), None, Some(300)).is_err());

        // The same with LCS contexts.
        let lcs = |rank: u32| rank % 7;
        let stored = chains.to_stored_with(Some(&lcs), max_direct_suffix_len(Some(300)), Some(300));
        assert_eq!(ChainCsr::from_stored_with(stored, Some(&lcs), Some(300)).unwrap(), chains);
        let stored = chains.to_stored_with(Some(&lcs), max_direct_suffix_len(Some(300)), Some(300));
        assert!(ChainCsr::from_stored_with(stored, Some(&lcs), None).is_err());
    }

    #[test]
    fn walks_are_bounded_by_the_contig_length() {
        // A corrupt graph: the source edge says 3 bases, then a node that loops to itself.
        let csr = Csr {
            node_weights: vec![SOURCE_COLEX, SINK_COLEX, 7],
            row_ptrs: vec![0, 1, 1, 2],
            column_indices: vec![2, 2],
            suffix_lens: vec![3, 1],
        };
        let choices = Choices { streams: vec![Default::default()] };
        assert!(csr.walk_path(&choices, 0).is_err());
        let chains = ChainCsr::from_csr(&csr);
        assert!(chains.walk_path(&choices, 0).is_err());
    }

    /// Random walks over many paths: chains with many back-references.
    fn many_paths() -> Vec<Vec<(usize, Range<usize>)>> {
        let mut x = 12345_u64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        (0..60).map(|_| (0..40).map(|_| { let r = (next() % 300) as usize; (1 + (next() % 5) as usize, r..r + 1) }).collect()).collect()
    }

    #[test]
    fn targets_roundtrip() {
        for dicts in [dicts(), many_paths()] {
            let (chains, choices) = builder(&dicts).to_chains();
            let decoded = ChainCsr::from_stored(chains.to_stored()).unwrap();
            assert_eq!(decoded, chains);
            for (color, d) in dicts.iter().enumerate() {
                assert_eq!(decoded.walk_path(&choices, color as u32).unwrap().1, expected(d));
            }
            // Corrupt targets give an error, not a panic.
            let mut bad = chains.to_stored();
            for b in bad.targets.iter_mut() { *b = !*b; }
            let _ = ChainCsr::from_stored(bad);
        }
    }

    #[test]
    fn stored_graph_roundtrip() {
        let dicts = dicts();
        for collapse in [false, true] {
            let (graph, choices) = builder(&dicts).collapse_chains(collapse).finish();
            let bytes = postcard::to_allocvec(&graph).unwrap();
            let graph: StoredGraph = postcard::from_bytes(&bytes).unwrap();
            assert_eq!(matches!(graph, StoredGraph::Chains(_)), collapse);
            for (color, d) in dicts.iter().enumerate() {
                let (_, decoded) = graph.walk_path(&choices, color as u32).unwrap();
                assert_eq!(decoded, expected(d));
            }
        }
    }
}
