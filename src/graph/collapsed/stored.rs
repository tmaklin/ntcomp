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
//! Stored form of [ChainCsr]: [StoredChains].
//!
//! The sections are the out-degrees, chain lengths and contig lengths
//! ([topology](super::topology)), the nodes and the chain targets
//! ([nodes](super::nodes)) and the suffix lengths
//! ([suffix_lens](super::suffix_lens)), each coded on its own.
//!
//! The nodes are coded in *node order*: the junctions in their breadth-first
//! numbering; for each junction its chains in order; along a chain, each
//! inner node, then the chain's target (with the index) or, without it, the
//! target junction if this chain is the first to reach it (see
//! [nodes](super::nodes)). The dummy source and sink are always junctions 0
//! and 1 and are not stored. The suffix lengths are coded in the same order,
//! interleaved: the length of an edge before the node it leads to.
//!
use crate::codecs::arith::{Decoder, Encoder};
use crate::codecs::bits::{BitReader, BitWriter};
use crate::graph::{GraphDecodeError, GraphWalk, SINK_COLEX, SINK_NODE, SOURCE_COLEX, SOURCE_NODE};
use crate::walk::{IndexWalk, WalkIndex};
use sbwt::StreamingIndex;

use super::nodes::{code_targets, NodeReader, NodeWriter, Prepared, TargetModels, WalkModels};
use rayon::prelude::*;
use super::suffix_lens::{max_direct_suffix_len, SuffixModel};
use super::topology::code_topology;
use super::ChainCsr;

/// Stored form of [ChainCsr], see the module docs.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct StoredChains {
    /// Number of junction nodes, including the dummy source and sink.
    pub n_junctions: u32,

    /// Every stored colex rank is below this (largest rank + 1).
    pub colex_bound: u32,

    /// Colex ranks of the stored (not derived) nodes in node order,
    /// truncated binary for `colex_bound`, see [nodes](super::nodes).
    pub node_ranks: Vec<u8>,

    /// Whether nodes are derived by walks in the index (the graph was stored
    /// with a [WalkIndex]); then decoding needs the index.
    pub derived: bool,

    /// Flags and branch choices of the derived nodes, arithmetic coded, see
    /// [nodes](super::nodes).
    pub node_walks: Vec<u8>,


    /// Largest LCS context of the suffix length models; 0: one context for
    /// all. Otherwise the context is min(LCS(source rank), lcs_cap).
    pub lcs_cap: u16,

    /// D: suffix lengths up to this are coded directly, see
    /// [suffix_lens](super::suffix_lens).
    pub max_direct: u32,

    /// Other suffix lengths, arithmetic coded, see
    /// [suffix_lens](super::suffix_lens).
    pub suffix_lens: Vec<u8>,

    /// Out-degree of every junction except the sink, the number of inner
    /// nodes of every chain and the contig lengths (suffix lengths of the
    /// edges out of the dummy source), arithmetic coded, see
    /// [topology](super::topology).
    pub topology: Vec<u8>,

    /// Target of every chain, arithmetic coded, see [nodes](super::nodes); empty
    /// when the nodes are derived (the targets are in `node_walks`).
    pub targets: Vec<u8>,
}

impl ChainCsr {
    /// Converts into the stored form.
    ///
    /// Panics if the junctions are not numbered in breadth-first order.
    pub fn to_stored(
        &self,
    ) -> StoredChains {
        self.to_stored_with(None, 0, None, false)
    }

    /// Same as [ChainCsr::to_stored], with the `index`: one suffix length
    /// code per LCS value of the edge's source rank (capped at `lcs_cap`),
    /// and if `derive`, nodes derived by walks in the index (which needs its
    /// LCS array); and direct suffix length symbols for _k_-mer size `k`. See
    /// the module docs, [suffix_lens](super::suffix_lens) and
    /// [nodes](super::nodes).
    pub fn to_stored_with(
        &self,
        index: Option<WalkIndex>,
        lcs_cap: u32,
        k: Option<u32>,
        derive: bool,
    ) -> StoredChains {
        let lcs_cap = if index.is_some() { lcs_cap.max(1) } else { 0 };
        let max_direct = max_direct_suffix_len(k);
        let n_junctions = self.junction_weights.len();

        let (sink_start, sink_end) = self.row(SINK_NODE as usize);
        assert_eq!(sink_start, sink_end, "sink has outgoing chains");
        let degrees: Vec<u64> = (0..n_junctions).map(|j| { let (start, end) = self.row(j); (end - start) as u64 }).collect();
        assert!(degrees.iter().enumerate().all(|(j, &d)| d > 0 || j == SINK_NODE as usize), "a junction has no outgoing chains");
        let inner: Vec<u64> = (0..self.n_chains()).map(|chain| self.inner_nodes(chain).len() as u64).collect();
        // Chains out of the source come first; their first edge holds the contig length.
        let n_source_chains = self.row(SOURCE_NODE as usize).1;
        let contig_lens: Vec<u64> = (0..n_source_chains).map(|chain| self.chain_suffix_lens(chain)[0] as u64).collect();
        let mut topology = Encoder::new();
        code_topology(&mut topology, n_junctions, Some((&degrees, &inner, &contig_lens))).expect("topology encode");
        let mut targets = Encoder::new();
        let (_, next_new) = code_targets(&mut targets, self.n_chains(), Some(&self.targets)).expect("targets encode");
        assert_eq!(next_new as usize, n_junctions, "junctions are not numbered in breadth-first order");

        // Real colex ranks only; the dummy source and sink are junctions 0 and 1.
        assert_eq!(self.junction_weights[..2], [SOURCE_COLEX, SINK_COLEX], "junctions 0 and 1 are not the source and sink");
        let real_junctions = &self.junction_weights[2..];
        let max_colex = real_junctions.iter().chain(&self.inner_weights).copied().max().unwrap_or(0);
        assert!(max_colex < SOURCE_COLEX, "colex rank collides with dummy node ranks");
        let colex_bound = max_colex + 1;

        let walk = derive.then(|| {
            let index = index.expect("deriving nodes needs the index");
            let streaming = StreamingIndex::new(index.sbwt, index.lcs.expect("deriving nodes needs the LCS array"));
            (index, streaming)
        });
        let derived = walk.is_some();
        let mut nodes = NodeWriter { walk, ranks: BitWriter::new(), walks: Encoder::new(), models: WalkModels::new(), bound: colex_bound };
        // With walks, the targets are coded where each chain ends.
        let mut end_targets = TargetModels::new();
        let mut model = SuffixModel::new(lcs_cap as usize + 1, max_direct);
        let mut suffix_lens = Encoder::new();
        // The nodes and targets in node order (predecessor, node; None for the
        // sink as a target), prepared in parallel when they are derived.
        let mut order: Vec<(Option<(u32, u32)>, Option<u32>)> = Vec::new();
        if derived {
            for junction in 0..n_junctions {
                let (start, end) = self.row(junction);
                for chain in start..end {
                    let inner = self.inner_nodes(chain);
                    let lens = self.chain_suffix_lens(chain);
                    for (i, &len) in lens.iter().enumerate() {
                        let source = if i == 0 { self.junction_weights[junction] } else { inner[i - 1] };
                        let from_source = junction == SOURCE_NODE as usize && i == 0;
                        let pred = (!from_source).then_some((source, len));
                        let x = if i + 1 < lens.len() {
                            Some(inner[i])
                        } else {
                            let target = self.targets[chain];
                            (target != SINK_NODE).then(|| self.junction_weights[target as usize])
                        };
                        order.push((pred, x));
                    }
                }
            }
        }
        let prepared: Vec<Prepared> = order.par_iter().map(|&(pred, x)| NodeWriter::prepare(&nodes.walk, pred, x)).collect();
        drop(order);
        let mut prepared = prepared.into_iter();

        let mut reached = vec![false; n_junctions];
        reached[SOURCE_NODE as usize] = true;
        reached[SINK_NODE as usize] = true;
        for junction in 0..n_junctions {
            let (start, end) = self.row(junction);
            for chain in start..end {
                let inner = self.inner_nodes(chain);
                let lens = self.chain_suffix_lens(chain);
                for (i, &len) in lens.iter().enumerate() {
                    let source = if i == 0 { self.junction_weights[junction] } else { inner[i - 1] };
                    let from_source = junction == SOURCE_NODE as usize && i == 0;
                    if !from_source {
                        let context = index.map_or(0, |index| index.lcs_at(source as usize).min(lcs_cap as usize));
                        model.code(&mut suffix_lens, context, len).expect("suffix length encode");
                    }
                    if i + 1 < lens.len() {
                        if derived {
                            nodes.node_prepared(prepared.next().unwrap(), inner[i], true);
                        } else {
                            nodes.node(inner[i]);
                        }
                    } else {
                        let target = self.targets[chain] as usize;
                        let is_new = !reached[target];
                        reached[target] = true;
                        if derived {
                            nodes.end_prepared(prepared.next().unwrap(), target as u32, self.junction_weights[target], is_new, &mut end_targets);
                        } else if is_new {
                            nodes.node(self.junction_weights[target]);
                        }
                    }
                }
            }
        }

        StoredChains {
            n_junctions: n_junctions.try_into().expect("too many junctions"),
            colex_bound,
            node_ranks: nodes.ranks.finish(),
            derived,
            node_walks: if derived { nodes.walks.finish() } else { Vec::new() },
            lcs_cap: lcs_cap as u16,
            max_direct,
            suffix_lens: suffix_lens.finish(),
            topology: topology.finish(),
            targets: if derived { Vec::new() } else { targets.finish() },
        }
    }

    /// Converts from the stored form.
    pub fn from_stored(
        stored: StoredChains,
    ) -> Result<Self, GraphDecodeError> {
        Self::from_stored_with(stored, None, None)
    }

    /// Same as [ChainCsr::from_stored]; `index` is needed if the chains were
    /// stored with it (LCS contexts or derived nodes; its LCS array is
    /// optional, the values are otherwise computed from the SBWT), and `k`
    /// must be the _k_-mer size they were stored with (None for _k_ ≤ 255).
    pub fn from_stored_with(
        stored: StoredChains,
        index: Option<WalkIndex>,
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

        let mut coded = Decoder::new(&stored.topology);
        let (degrees, inner, contig_lens) = code_topology(&mut coded, n_junctions, None).ok_or_else(|| truncated("out-degrees, chain lengths and contig lengths"))?;
        if coded.overrun() {
            return Err(truncated("out-degrees, chain lengths and contig lengths"))
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

        let lcs_cap = stored.lcs_cap as u32;
        if index.is_none() && lcs_cap > 0 {
            return Err(GraphDecodeError("stored chains: the suffix lengths need the index".to_string()))
        }
        if index.is_none() && stored.derived {
            return Err(GraphDecodeError("stored chains: the derived nodes need the index".to_string()))
        }
        let lcs = if lcs_cap > 0 { index } else { None };
        let walk = if stored.derived { index } else { None };

        // The targets set the node order: stored first, or with walks coded
        // where each chain ends.
        let mut targets = if stored.derived {
            vec![0_u32; n_chains]
        } else {
            let mut coded = Decoder::new(&stored.targets);
            let (targets, next_new) = code_targets(&mut coded, n_chains, None)?;
            if next_new as usize != n_junctions || coded.overrun() {
                return Err(truncated("targets"))
            }
            targets
        };

        // Nodes and suffix lengths in node order, see the module docs.
        let mut nodes = NodeReader {
            walk,
            ranks: BitReader::new(&stored.node_ranks),
            walks: Decoder::new(&stored.node_walks),
            models: WalkModels::new(),
            bound: colex_bound,
        };
        let mut end_targets = TargetModels::new();
        // Colex rank -> junction, for the junctions reached so far.
        let mut junction_of: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
        let mut junction_weights = vec![0_u32; n_junctions];
        junction_weights[SOURCE_NODE as usize] = SOURCE_COLEX;
        junction_weights[SINK_NODE as usize] = SINK_COLEX;
        let mut reached = vec![false; n_junctions];
        reached[SOURCE_NODE as usize] = true;
        reached[SINK_NODE as usize] = true;
        let mut inner_weights = vec![0_u32; n_inner];
        let mut model = SuffixModel::new(lcs_cap as usize + 1, max_direct);
        let mut coded = Decoder::new(&stored.suffix_lens);
        let mut suffix_lens = Vec::with_capacity(n_inner + n_chains);
        for junction in 0..n_junctions {
            if !reached[junction] {
                return Err(truncated("targets"))
            }
            for chain in row_ptrs[junction] as usize..row_ptrs[junction + 1] as usize {
                let start = chain_ptrs[chain] as usize;
                let n_edges = chain_ptrs[chain + 1] as usize - start + 1;
                for i in 0..n_edges {
                    let source = if i == 0 { junction_weights[junction] } else { inner_weights[start + i - 1] };
                    let from_source = junction == SOURCE_NODE as usize && i == 0;
                    let len = if from_source {
                        contig_lens[chain] as u32
                    } else {
                        let context = lcs.map_or(0, |index| index.lcs_at(source as usize).min(lcs_cap as usize));
                        model.code(&mut coded, context, 0).ok_or_else(|| truncated("suffix lengths"))?
                    };
                    suffix_lens.push(len);
                    let pred = (!from_source).then_some((source, len));
                    if i + 1 < n_edges {
                        inner_weights[start + i] = nodes.node(pred, true).ok_or_else(|| truncated("nodes"))?;
                    } else if stored.derived {
                        let (target, rank) = nodes.end(pred, &mut end_targets, &junction_of).ok_or_else(|| truncated("targets"))?;
                        if target as usize >= n_junctions {
                            return Err(truncated("targets"))
                        }
                        targets[chain] = target;
                        if let Some(rank) = rank {
                            reached[target as usize] = true;
                            junction_weights[target as usize] = rank;
                            junction_of.insert(rank, target);
                        }
                    } else {
                        let target = targets[chain] as usize;
                        if !reached[target] {
                            reached[target] = true;
                            junction_weights[target] = nodes.node(pred, false).ok_or_else(|| truncated("nodes"))?;
                        }
                    }
                }
            }
        }
        if coded.overrun() {
            return Err(truncated("suffix lengths"))
        }
        if (stored.derived && (nodes.walks.overrun() || end_targets.next_new as usize != n_junctions))
            || nodes.ranks.bits_read().div_ceil(8) != stored.node_ranks.len() {
            return Err(truncated("nodes"))
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
}

#[cfg(test)]
mod tests {
    use crate::graph::*;
    use super::super::tests::{builder, dicts, expected, many_paths};
    use super::super::{max_direct_suffix_len, ChainCsr, MIN_DIRECT_SUFFIX_LEN};
    use crate::nlz::Entry;
    use crate::testing::{dictionary, dictionary_with_flags};
    use crate::testing::related_index;

    use core::ops::Range;

    #[test]
    fn stored_chains_roundtrip() {
        let dicts = dicts();
        let (chains, _) = builder(&dicts).finish();
        assert_eq!(chains.junction_weights[..2], [SOURCE_COLEX, SINK_COLEX]);

        let stored = chains.to_stored();
        assert_eq!(ChainCsr::from_stored(stored).unwrap(), chains);

        for field in 0..9 {
            let mut stored = chains.to_stored();
            match field {
                0 => stored.topology.clear(),
                1 => { stored.topology.pop(); },
                2 => stored.suffix_lens.clear(),
                3 => { stored.topology.truncate(stored.topology.len() / 2); },
                4 => { stored.node_ranks.pop(); },
                5 => { stored.node_ranks.push(0); },
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
        let (chains, _) = builder(&dicts()).finish();
        let stored = chains.to_stored();
        assert_eq!(stored.colex_bound, 96);
        assert_eq!(stored.n_junctions as usize, chains.junction_weights.len());
        let bits = |ranks: &[u32]| ranks.iter().map(|&r| if r < 32 { 6 } else { 7 }).sum::<usize>().div_ceil(8);
        // Without the index every node is stored.
        let all: Vec<u32> = chains.junction_weights[2..].iter().chain(&chains.inner_weights).copied().collect();
        assert_eq!(stored.node_ranks.len(), bits(&all));
        assert!(!stored.derived && stored.node_walks.is_empty());

        // Ranks of 0 only (no bits), a power-of-two bound, and ranks just below the dummy ranks.
        for (ranks, bound) in [(vec![0, 0], 1), (vec![5, 63, 0], 64), (vec![1, SOURCE_COLEX - 1, 7], SOURCE_COLEX)] {
            let d: Vec<Vec<(usize, Range<usize>)>> = vec![ranks.iter().map(|&r| (3, r as usize..r as usize + 1)).collect()];
            let (chains, choices) = builder(&d).finish();
            let stored = chains.to_stored();
            assert_eq!(stored.colex_bound, bound);
            if bound == 1 {
                assert!(stored.node_ranks.is_empty());
            }
            let decoded = ChainCsr::from_stored(stored).unwrap();
            assert_eq!(decoded, chains);
            let path = decoded.walk_path(&choices, 0).unwrap();
            assert_eq!(path.entries.iter().map(|e| e.colex_rank).collect::<Vec<u32>>(), ranks);
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
            builder.add_dictionary(&dictionary_with_flags(d, f));
        }
        let (chains, choices) = builder.finish();

        let bytes = postcard::to_allocvec(&chains).unwrap();
        let decoded: ChainCsr = postcard::from_bytes(&bytes).unwrap();
        assert_eq!(decoded, chains);

        for (color, (d, implicit)) in dicts.iter().zip(&expected_implicit).enumerate() {
            let path = decoded.walk_path(&choices, color as u32).unwrap();
            assert_eq!(path.sequence_len, d.iter().map(|x| x.0).sum::<usize>());
            // Implicit entries come back with their length unknown.
            let expected: Vec<Entry> = d.iter().zip(implicit).map(|((s, r), &implicit)| {
                Entry { colex_rank: r.start as u32, len: if implicit { 0 } else { *s as u32 }, implicit }
            }).collect();
            assert_eq!(path.entries, expected);
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
            builder.add_dictionary(&dictionary(d));
        }
        let (chains, _) = builder.finish();

        // For k <= 255 the coding is the same as without k.
        for k in [Some(31), Some(255)] {
            assert_eq!(postcard::to_allocvec(&chains.to_stored_with(None, 0, k, false)).unwrap(), postcard::to_allocvec(&chains.to_stored()).unwrap());
        }
        assert_eq!(max_direct_suffix_len(Some(63)), MIN_DIRECT_SUFFIX_LEN);
        assert_eq!(max_direct_suffix_len(Some(300)), 300);

        // With k = 300 every length up to 300 is coded directly, and the chains
        // only decode with the same k.
        let stored = chains.to_stored_with(None, 0, Some(300), false);
        assert_eq!(stored.max_direct, 300);
        assert_eq!(chains.to_stored().max_direct, MIN_DIRECT_SUFFIX_LEN);
        assert_eq!(ChainCsr::from_stored_with(chains.to_stored_with(None, 0, Some(300), false), None, Some(300)).unwrap(), chains);
        assert!(ChainCsr::from_stored_with(chains.to_stored_with(None, 0, Some(300), false), None, None).is_err());
        assert!(ChainCsr::from_stored_with(chains.to_stored(), None, Some(300)).is_err());
    }

    #[test]
    fn targets_roundtrip() {
        for dicts in [dicts(), many_paths()] {
            let (chains, choices) = builder(&dicts).finish();
            let decoded = ChainCsr::from_stored(chains.to_stored()).unwrap();
            assert_eq!(decoded, chains);
            for (color, d) in dicts.iter().enumerate() {
                assert_eq!(decoded.walk_path(&choices, color as u32).unwrap(), expected(d));
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
        let (graph, choices) = builder(&dicts).finish();
        let mut bytes = graph.to_bytes(None, None).unwrap();
        let decoded = ChainCsr::from_bytes(&bytes, None, None).unwrap();
        assert_eq!(decoded, graph);
        for (color, d) in dicts.iter().enumerate() {
            assert_eq!(decoded.walk_path(&choices, color as u32).unwrap(), expected(d));
        }
        // The plain layout's tag and an empty section are errors.
        bytes[0] = 0;
        assert!(ChainCsr::from_bytes(&bytes, None, None).is_err());
        assert!(ChainCsr::from_bytes(&[], None, None).is_err());
    }

    #[test]
    fn lcs_context_suffix_codes_roundtrip() {
        let (seqs, index) = related_index(15);
        let (sbwt, lcs) = (index.sbwt(), index.lcs());
        let mut builder = crate::graph::GraphBuilder::new();
        for seq in &seqs {
            builder.add_dictionary(&crate::nlz::parse(seq, &index).unwrap());
        }
        let (chains, _) = builder.finish();
        let from_array = crate::walk::WalkIndex { sbwt, lcs: Some(lcs) };
        let from_sbwt = crate::walk::WalkIndex { sbwt, lcs: None };
        for cap in [1, 8, 255] {
            let stored = chains.to_stored_with(Some(from_array), cap, Some(15), false);
            let bytes = postcard::to_allocvec(&stored).unwrap();
            // Encoded with the array, decoded with values computed from the SBWT (and vice versa).
            let decoded = crate::graph::collapsed::ChainCsr::from_stored_with(postcard::from_bytes(&bytes).unwrap(), Some(from_sbwt), Some(15)).unwrap();
            assert_eq!(decoded, chains, "cap {}", cap);
            let stored = chains.to_stored_with(Some(from_sbwt), cap, Some(15), false);
            assert_eq!(crate::graph::collapsed::ChainCsr::from_stored_with(stored, Some(from_array), Some(15)).unwrap(), chains);
            // The index is required to decode.
            let stored: crate::graph::collapsed::StoredChains = postcard::from_bytes(&bytes).unwrap();
            assert!(crate::graph::collapsed::ChainCsr::from_stored_with(stored, None, Some(15)).is_err());
        }
        // Without contexts nothing changes.
        assert_eq!(crate::graph::collapsed::ChainCsr::from_stored(chains.to_stored()).unwrap(), chains);
    }

    #[test]
    fn nodes_derived_by_walks_roundtrip() {
        use crate::graph::GraphWalk;
        let (seqs, index) = related_index(15);
        let (sbwt, lcs) = (index.sbwt(), index.lcs());
        let mut builder = crate::graph::GraphBuilder::new();
        for seq in &seqs {
            builder.add_dictionary(&crate::nlz::parse(seq, &index).unwrap());
        }
        let (chains, _) = builder.finish();
        let walk = crate::graph::collapsed::WalkIndex { sbwt, lcs: Some(lcs) };
        let plain = chains.to_stored_with(Some(walk), 255, Some(15), false);
        let stored = chains.to_stored_with(Some(walk), 255, Some(15), true);
        assert!(stored.derived && !plain.derived);
        // Some nodes are derived: fewer ranks are stored.
        assert!(stored.node_ranks.len() < plain.node_ranks.len(), "{} vs {}", stored.node_ranks.len(), plain.node_ranks.len());
        // With walks the targets are coded where each chain ends, some of
        // them junctions reached before.
        assert!(stored.targets.is_empty() && !plain.targets.is_empty());
        let mut reached = vec![false; chains.junction_weights.len()];
        let mut reached_before = 0;
        for junction in 0..chains.junction_weights.len() {
            let (start, end) = chains.row(junction);
            for chain in start..end {
                let target = chains.targets[chain] as usize;
                reached_before += reached[target] as usize;
                reached[target] = true;
            }
        }
        assert!(reached_before > 0, "no chain ends at a junction reached before");
        let bytes = postcard::to_allocvec(&stored).unwrap();
        // Decoded with the LCS array, and with LCS values computed from the SBWT (no array).
        let decoded = crate::graph::collapsed::ChainCsr::from_stored_with(postcard::from_bytes(&bytes).unwrap(), Some(walk), Some(15)).unwrap();
        assert_eq!(decoded, chains);
        let no_array = crate::graph::collapsed::WalkIndex { sbwt, lcs: None };
        let decoded = crate::graph::collapsed::ChainCsr::from_stored_with(postcard::from_bytes(&bytes).unwrap(), Some(no_array), Some(15)).unwrap();
        assert_eq!(decoded, chains);
        // The derived nodes need the index.
        assert!(crate::graph::collapsed::ChainCsr::from_stored_with(postcard::from_bytes(&bytes).unwrap(), None, Some(15)).is_err());
        // A corrupt walk stream is an error, not a panic.
        let mut bad: crate::graph::collapsed::StoredChains = postcard::from_bytes(&bytes).unwrap();
        bad.node_walks.truncate(bad.node_walks.len() / 2);
        let _ = crate::graph::collapsed::ChainCsr::from_stored_with(bad, Some(walk), Some(15));
    }
}
