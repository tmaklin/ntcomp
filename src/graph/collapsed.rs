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
//! A node with exactly one incoming and one outgoing edge is an *inner* node;
//! the dummy source and sink and every other node are *junctions*. Following
//! edges from a junction through inner nodes always ends at a junction, so the
//! graph can be stored as chains between junctions.
//!
//! A chain stores the colex ranks of its inner nodes and the suffix lengths of
//! its edges, and its colour set only once. An inner node passes every
//! colour that enters it on to its only outgoing edge, so all edges in a chain
//! have the same colours.
//!
//! Colours only have a choice at junctions, so the choice stream has the same
//! number of entries as for the uncollapsed graph. The stored ranks refer to the
//! order of the chains in the junction's row.
//!
use super::{
    Csr,
    Choices,
    GraphDecodeError,
    SINK_NODE,
    SOURCE_NODE,
    take_edge,
};

/// Graph of chains between junctions, in compressed sparse row format.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct ChainCsr {
    /// Colex ranks of the junctions; junction 0 is the source and 1 the sink.
    pub junction_weights: Vec<u32>,

    /// Chains leaving junction `j` are `row_ptrs[j]..row_ptrs[j + 1]`.
    pub row_ptrs: Vec<u32>,

    /// Junction that each chain ends at.
    pub targets: Vec<u32>,

    /// Inner nodes of chain `c` are `inner_weights[chain_ptrs[c]..chain_ptrs[c + 1]]`.
    pub chain_ptrs: Vec<u32>,

    /// Colex ranks of the inner nodes, concatenated over all chains.
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

        // Node index -> junction index, u32::MAX for inner nodes. The source
        // and sink are nodes 0 and 1, so they become junctions 0 and 1.
        let mut junction_index = vec![u32::MAX; n_nodes];
        let mut junction_weights = Vec::new();
        for node in 0..n_nodes {
            let (start, end) = csr.row(node);
            let is_dummy = node == SOURCE_NODE as usize || node == SINK_NODE as usize;
            if is_dummy || in_degree[node] != 1 || end - start != 1 {
                junction_index[node] = junction_weights.len() as u32;
                junction_weights.push(csr.node_weights[node]);
            }
        }

        let mut row_ptrs = Vec::with_capacity(junction_weights.len() + 1);
        let mut targets = Vec::new();
        let mut chain_ptrs = vec![0_u32];
        let mut inner_weights = Vec::new();
        let mut suffix_lens = Vec::with_capacity(n_edges);

        row_ptrs.push(0);
        for node in 0..n_nodes {
            if junction_index[node] == u32::MAX {
                continue;
            }
            let (start, end) = csr.row(node);
            for pos in start..end {
                suffix_lens.push(csr.suffix_lens[pos]);

                let mut next = csr.column_indices[pos] as usize;
                while junction_index[next] == u32::MAX {
                    let edge = csr.row(next).0;
                    inner_weights.push(csr.node_weights[next]);
                    suffix_lens.push(csr.suffix_lens[edge]);
                    next = csr.column_indices[edge] as usize;
                }

                targets.push(junction_index[next]);
                chain_ptrs.push(inner_weights.len() as u32);
            }
            row_ptrs.push(targets.len() as u32);
        }

        ChainCsr {
            junction_weights,
            row_ptrs,
            targets,
            chain_ptrs,
            inner_weights,
            suffix_lens,
        }
    }

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

    /// Walks the path of `color` from the source to the sink.
    ///
    /// Same as [super::walk_path] for the uncollapsed graph; `max_steps` bounds
    /// the number of uncollapsed edges.
    pub fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
        max_steps: usize,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError> {
        let mut reader = choices
            .reader(color)
            .ok_or_else(|| GraphDecodeError(format!("no choice stream for colour {}", color)))?;
        let mut junction = SOURCE_NODE as usize;
        let mut sequence_length = 0_u32;
        let mut dictionary: Vec<(u32, u32)> = Vec::new();
        let mut steps = 0;

        loop {
            let chain = take_edge(self.row(junction), &mut reader, color, junction)?;
            let inner = self.inner_nodes(chain);
            let suffix_lens = self.chain_suffix_lens(chain);

            steps += suffix_lens.len();
            if steps > max_steps {
                return Err(GraphDecodeError(format!("colour {} did not reach the sink in {} steps", color, max_steps)))
            }

            if junction == SOURCE_NODE as usize {
                sequence_length = suffix_lens[0];
            } else {
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
                for i in 0..self.inner_nodes(chain).len() {
                    let node = n_junctions + self.chain_ptrs[chain] + i as u32;
                    edges.push((prev, node, suffix_lens[i]));
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
