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
//! Building the graph from the contigs' NLZ dictionaries ([nlz::parse](crate::nlz::parse)),
//! one path per colour. Choosing the strands and editing the contigs before
//! parsing them is the archive's ([crate::archive::ArchiveWriter]).
//!
use core::ops::Range;

use indexmap::IndexSet;

use super::choices::{ChoiceWriter, Choices};
use super::collapsed::ChainCsr;
use crate::nlz::IMPLICIT_SUFFIX_LEN;
use super::{Csr, GraphWalk, StoredGraph, SINK_COLEX, SINK_NODE, SOURCE_COLEX, SOURCE_NODE};

/// Builds the graph and remembers the path each colour takes.
pub struct GraphBuilder {

    /// colex rank -> node index (position in the set)
    nodes: IndexSet<u32>,

    /// (source node, target node, suffix length); edge id is the position in the set
    edges: IndexSet<(u32, u32, u32)>,

    /// Per edge id: true if the suffix length is implicit, see [SuffixLen](crate::nlz::SuffixLen)
    implicit: Vec<bool>,

    /// Per colour: edge ids (indexes into `edges`) in traversal order
    paths: Vec<Vec<usize>>,

    /// Store linear chains as single edges in [GraphBuilder::finish]
    collapse_chains: bool,
}

impl Default for GraphBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl GraphBuilder {
    pub fn new() -> Self {
        let mut nodes = IndexSet::new();
        nodes.insert(SOURCE_COLEX);
        nodes.insert(SINK_COLEX);

        GraphBuilder {
            nodes,
            edges: IndexSet::new(),
            implicit: Vec::new(),
            paths: Vec::new(),
            collapse_chains: false,
        }
    }

    /// Collapse linear chains into single edges in [GraphBuilder::finish].
    pub fn collapse_chains(
        mut self,
        collapse: bool,
    ) -> Self {
        self.collapse_chains = collapse;
        self
    }

    fn insert_node(
        &mut self,
        colex_rank: usize
    ) -> u32 {
        let colex_rank: u32 = colex_rank.try_into().expect("colex rank does not fit in u32");
        assert!(colex_rank < SOURCE_COLEX, "colex rank collides with dummy node ranks");
        self.nodes.insert_full(colex_rank).0 as u32
    }

    fn edge(
        &mut self,
        source: u32,
        target: u32,
        weight: u32,
        implicit: bool,
    ) -> usize {
        let (id, inserted) = self.edges.insert_full((source, target, weight));
        if inserted {
            self.implicit.push(false);
        }
        // Implicit if any colour saw the extension stop in the SBWT: the length
        // then only depends on the _k_-mer, not on where the colour uses it.
        self.implicit[id] |= implicit;
        id
    }

    /// Adds the dictionary of the next colour.
    ///
    /// `dictionary` is the output of [nlz::parse](crate::nlz::parse). The edge out of
    /// the source node stores `sequence_length`.
    pub fn add_path(
        &mut self,
        sequence_length: usize,
        dictionary: &[(usize, Range<usize>)],
    ) -> u32 {
        self.add_path_with_flags(sequence_length, dictionary, &vec![false; dictionary.len()])
    }

    /// Same as [GraphBuilder::add_path] but stores the suffix length of the
    /// entries flagged in `implicit` as [SuffixLen::Implicit](crate::nlz::SuffixLen::Implicit).
    ///
    /// `dictionary` and `implicit` are the output of [nlz::parse](crate::nlz::parse).
    pub fn add_path_with_flags(
        &mut self,
        sequence_length: usize,
        dictionary: &[(usize, Range<usize>)],
        implicit: &[bool],
    ) -> u32 {
        assert_eq!(dictionary.len(), implicit.len());
        let color: u32 = self.paths.len().try_into().unwrap();
        let mut path = Vec::with_capacity(dictionary.len() + 1);

        let mut prev = SOURCE_NODE;
        let mut prev_weight: u32 = sequence_length.try_into().expect("sequence length does not fit in u32");
        let mut prev_implicit = false;
        for ((suffix_len, colex), &is_implicit) in dictionary.iter().zip(implicit) {
            let node = self.insert_node(colex.start);
            path.push(self.edge(prev, node, prev_weight, prev_implicit));
            prev = node;
            prev_weight = (*suffix_len).try_into().unwrap();
            assert!(prev_weight != 0, "empty dictionary entry");
            prev_implicit = is_implicit;
        }
        path.push(self.edge(prev, SINK_NODE, prev_weight, prev_implicit));

        self.paths.push(path);
        color
    }

    /// Converts into the stored representation (graph + choice stream),
    /// collapsing linear chains if [GraphBuilder::collapse_chains] was set.
    pub fn finish(self) -> (StoredGraph, Choices) {
        if self.collapse_chains {
            let (chains, choices) = self.to_chains();
            (StoredGraph::Chains(chains), choices)
        } else {
            let (csr, choices) = self.to_csr();
            (StoredGraph::Csr(csr), choices)
        }
    }

    /// Converts into the stored representation (graph + choice stream).
    pub fn to_csr(self) -> (Csr<u32>, Choices) {
        let (csr, csr_pos) = self.build_csr();
        let choices = self.choices(&csr, &csr_pos);
        (csr, choices)
    }

    /// Converts into the stored representation with linear chains collapsed.
    ///
    /// Only junctions have more than one outgoing edge and their chains keep
    /// the CSR order of their first edges, so the choice stream is the same as
    /// for [GraphBuilder::to_csr].
    pub fn to_chains(self) -> (ChainCsr, Choices) {
        let (csr, csr_pos) = self.build_csr();
        let choices = self.choices(&csr, &csr_pos);
        (ChainCsr::from_csr(&csr), choices)
    }

    /// Records which edge each colour takes at every node with more than one outgoing edge.
    fn choices(
        &self,
        csr: &Csr<u32>,
        csr_pos: &[u32],
    ) -> Choices {
        let streams = self.paths.iter().map(|path| {
            let mut writer = ChoiceWriter::new();
            for &id in path {
                let source = self.edges.get_index(id).unwrap().0 as usize;
                let (start, end) = csr.row(source);
                if end - start > 1 {
                    writer.push((csr_pos[id] as usize - start) as u32, (end - start) as u32);
                }
            }
            writer.finish()
        }).collect();
        Choices { streams }
    }

    /// Builds the CSR and the CSR position of every edge id.
    fn build_csr(&self) -> (Csr<u32>, Vec<u32>) {
        let n_nodes = self.nodes.len();
        let n_edges = self.edges.len();

        // Number of times the colours take each edge.
        let mut uses = vec![0_u32; n_edges];
        for path in &self.paths {
            for &id in path {
                uses[id] += 1;
            }
        }

        // Edge order: by source, then the most used edges first, then target
        // and suffix length. The choice streams store ranks in this order, and
        // rank 0 being the most common choice makes them compress well.
        let mut order: Vec<usize> = (0..n_edges).collect();
        order.sort_unstable_by_key(|&id| {
            let &(source, target, weight) = self.edges.get_index(id).unwrap();
            (source, std::cmp::Reverse(uses[id]), target, weight)
        });

        let mut csr_pos = vec![0_u32; n_edges];
        let mut row_ptrs = vec![0_u32; n_nodes + 1];
        let mut column_indices = Vec::with_capacity(n_edges);
        let mut suffix_lens = Vec::with_capacity(n_edges);

        for (pos, &id) in order.iter().enumerate() {
            let &(source, target, weight) = self.edges.get_index(id).unwrap();
            csr_pos[id] = pos as u32;
            row_ptrs[source as usize + 1] += 1;
            column_indices.push(target);
            suffix_lens.push(if self.implicit[id] { IMPLICIT_SUFFIX_LEN } else { weight });
        }
        for i in 0..n_nodes {
            row_ptrs[i + 1] += row_ptrs[i];
        }

        let csr = Csr {
            node_weights: self.nodes.iter().copied().collect(),
            row_ptrs,
            column_indices,
            suffix_lens,
        };

        (csr, csr_pos)
    }
}

