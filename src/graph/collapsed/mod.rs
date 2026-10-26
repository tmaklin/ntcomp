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
//!
//! ## Stored form
//!
//! [StoredChains] ([stored]) holds the graph in sections, each coded on its
//! own:
//! - the out-degrees of the junctions, the chain lengths and the contig
//!   lengths ([topology]);
//! - the nodes, as colex ranks or derived by walks in the index, and the
//!   chain targets ([nodes]);
//! - the suffix lengths ([suffix_lens]).
//!
use crate::nlz::{Dictionary, Entry};
use super::{
    stored_entry,
    Csr,
    Choices,
    GraphDecodeError,
    SINK_NODE,
    SOURCE_NODE,
    take_edge,
    GraphWalk,
};

pub mod nodes;
pub mod stored;
pub mod suffix_lens;
pub mod topology;

pub use crate::walk::WalkIndex;
pub use nodes::TARGET_MODELLED_BITS;
pub use stored::StoredChains;
pub use suffix_lens::{max_direct_suffix_len, MIN_DIRECT_SUFFIX_LEN};

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


    /// Number of linear chains.
    pub fn n_chains(&self) -> usize {
        self.targets.len()
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

impl ChainCsr {
    /// The stored length of `color`'s contig: the suffix length of the
    /// chain its path takes out of the source.
    pub fn sequence_len(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<u64, GraphDecodeError> {
        let mut reader = choices
            .reader(color)
            .ok_or_else(|| GraphDecodeError(format!("no choice stream for colour {}", color)))?;
        let chain = take_edge(self.row(SOURCE_NODE as usize), &mut reader, color, SOURCE_NODE as usize)?;
        Ok(self.chain_suffix_lens(chain)[0] as u64)
    }
}

impl GraphWalk for ChainCsr {
    /// Number of junctions.
    fn n_nodes(&self) -> usize {
        self.row_ptrs.len().saturating_sub(1)
    }

    /// Number of chains.
    fn n_edges(&self) -> usize {
        self.n_chains()
    }

    /// Range of chains leaving junction `node`.
    fn row(&self, node: usize) -> (usize, usize) {
        (self.row_ptrs[node] as usize, self.row_ptrs[node + 1] as usize)
    }

    /// The junction `edge` (a chain) ends at.
    fn target(&self, edge: usize) -> usize {
        self.targets[edge] as usize
    }

    fn max_steps(&self, edge: usize) -> usize {
        (self.chain_suffix_lens(edge)[0] as usize).saturating_add(1)
    }

    /// Walks the path of `color` from the source to the sink.
    ///
    /// Same as the walk of the uncollapsed graph; the bound from [GraphWalk]
    /// counts uncollapsed edges.
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<Dictionary, GraphDecodeError> {
        let mut reader = choices
            .reader(color)
            .ok_or_else(|| GraphDecodeError(format!("no choice stream for colour {}", color)))?;
        let mut junction = SOURCE_NODE as usize;
        let mut sequence_length = 0_u32;
        let mut entries: Vec<Entry> = Vec::new();
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
                entries.push(stored_entry(self.junction_weights[junction], suffix_lens[0]));
            }
            entries.extend(inner.iter().zip(&suffix_lens[1..]).map(|(&colex_rank, &len)| stored_entry(colex_rank, len)));

            junction = self.targets[chain] as usize;
            if junction == SINK_NODE as usize {
                if !reader.is_finished() {
                    return Err(GraphDecodeError(format!("unused choices left for colour {}", color)))
                }
                return Ok(Dictionary { sequence_len: sequence_length as usize, entries })
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

    use core::ops::Range;
    use crate::testing::dictionary;

    pub(super) fn expected(dict: &[(usize, Range<usize>)]) -> crate::nlz::Dictionary {
        dictionary(dict)
    }

    // Shared linear stretches (50 -> 60 -> 70) and a repeat at node 10.
    pub(super) fn dicts() -> Vec<Vec<(usize, Range<usize>)>> {
        vec![
            vec![(5, 10..11), (1, 50..51), (2, 60..61), (3, 70..71), (5, 10..11), (4, 30..31), (5, 10..11), (3, 20..21)],
            vec![(5, 10..11), (4, 30..31), (5, 10..11), (1, 50..51), (2, 60..61), (3, 70..71), (2, 40..41)],
            vec![(6, 80..81), (1, 50..51), (2, 60..61), (3, 70..71), (2, 40..41), (7, 90..91)],
            vec![(9, 95..96)],
        ]
    }

    pub(super) fn builder(dicts: &[Vec<(usize, Range<usize>)>]) -> GraphBuilder {
        let mut builder = GraphBuilder::new();
        for d in dicts {
            builder.add_dictionary(&dictionary(d));
        }
        builder
    }

    #[test]
    fn chain_walk_recovers_paths() {
        let dicts = dicts();
        let (chains, choices) = builder(&dicts).finish();
        let csr = chains.to_csr();

        assert!(chains.n_chains() < csr.column_indices.len());
        assert_eq!(chains.suffix_lens.len(), csr.suffix_lens.len());

        let bytes = postcard::to_allocvec(&chains).unwrap();
        let chains: ChainCsr = postcard::from_bytes(&bytes).unwrap();

        for (color, d) in dicts.iter().enumerate().rev() {
            let decoded = chains.walk_path(&choices, color as u32).unwrap();
            assert_eq!(decoded, expected(d));
            assert_eq!(decoded.sequence_len, d.iter().map(|x| x.0).sum::<usize>());
        }
    }

    #[test]
    fn expanded_chains_decode_with_the_same_choices() {
        let dicts = dicts();
        let (chains, choices) = builder(&dicts).finish();
        let expanded = chains.to_csr();

        for (color, d) in dicts.iter().enumerate() {
            let decoded = expanded.walk_path(&choices, color as u32).unwrap();
            assert_eq!(decoded, expected(d));
        }
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

        let (chains, choices) = builder(&dicts).finish();
        let csr = chains.to_csr();
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
    pub(super) fn many_paths() -> Vec<Vec<(usize, Range<usize>)>> {
        let mut x = 12345_u64;
        let mut next = || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        (0..60).map(|_| (0..40).map(|_| { let r = (next() % 300) as usize; (1 + (next() % 5) as usize, r..r + 1) }).collect()).collect()
    }

}
