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
//! The graph in compressed sparse row format.
//!
#[cfg(feature = "graph-view")]
use petgraph::{
    EdgeType,
    graph::{
        Graph,
        NodeIndex,
    },
};

#[cfg(feature = "graph-view")]
use super::ColexGraphEdge;
use super::choices::Choices;
use super::paths::take_edge;
use super::{GraphDecodeError, GraphWalk, SINK_NODE, SOURCE_NODE};

/// Graph in compressed sparse row format; this is what gets stored.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct Csr<N> {
    pub node_weights: Vec<N>,
    pub row_ptrs: Vec<u32>,
    pub column_indices: Vec<u32>,
    pub suffix_lens: Vec<u32>,
}

impl<N: Copy> Csr<N> where u32: From<N> {
    /// Walks the path of `color` from the source to the sink, calling
    /// `visit(node, edge position)` for every edge taken.
    ///
    /// The walk is bounded by the contig length, see [GraphWalk].
    pub fn walk(
        &self,
        choices: &Choices,
        color: u32,
        mut visit: impl FnMut(usize, usize),
    ) -> Result<(), GraphDecodeError> {
        let mut reader = choices
            .reader(color)
            .ok_or_else(|| GraphDecodeError(format!("no choice stream for colour {}", color)))?;
        let mut node = SOURCE_NODE as usize;

        // Only the edge out of the source until its contig length is known.
        let mut max_edges = 1_usize;
        let mut n_edges = 0_usize;
        while n_edges < max_edges {
            let pos = take_edge(self.row(node), &mut reader, color, node)?;
            if node == SOURCE_NODE as usize {
                max_edges = (self.suffix_lens[pos] as usize).saturating_add(1);
            }
            n_edges += 1;
            visit(node, pos);

            node = self.column_indices[pos] as usize;
            if node == SINK_NODE as usize {
                if !reader.is_finished() {
                    return Err(GraphDecodeError(format!("unused choices left for colour {}", color)))
                }
                return Ok(())
            }
        }

        Err(GraphDecodeError(format!("colour {} did not reach the sink in {} edges", color, max_edges)))
    }

    /// Sorted colours using each edge, recovered by walking colours `0..n_colors`.
    pub fn colorsets(
        &self,
        choices: &Choices,
        n_colors: u32,
    ) -> Result<Vec<Vec<u32>>, GraphDecodeError> {
        let mut colorsets: Vec<Vec<u32>> = vec![Vec::new(); self.column_indices.len()];
        for color in 0..n_colors {
            self.walk(choices, color, |_, pos| {
                if colorsets[pos].last() != Some(&color) {
                    colorsets[pos].push(color);
                }
            })?;
        }
        Ok(colorsets)
    }

    /// Converts to a petgraph graph with `colorsets` from [Csr::colorsets] as the edge colours.
    #[cfg(feature = "graph-view")]
    pub fn to_petgraph<Ty: EdgeType>(
        self,
        colorsets: &[Vec<u32>],
    ) -> Graph<N, ColexGraphEdge, Ty>
    where
        N: Copy,
    {
        let mut graph = Graph::with_capacity(self.node_weights.len(), self.column_indices.len());

        for weight in &self.node_weights {
            graph.add_node(*weight);
        }

        for source_idx in 0..self.node_weights.len() {
            let (start, end) = self.row(source_idx);
            for (pos, colorset) in colorsets.iter().enumerate().take(end).skip(start) {
                graph.add_edge(
                    NodeIndex::new(source_idx),
                    NodeIndex::new(self.column_indices[pos] as usize),
                    ColexGraphEdge {
                        weight: self.suffix_lens[pos],
                        colors: colorset.clone(),
                    },
                );
            }
        }

        graph
    }
}

impl<N: Copy> GraphWalk for Csr<N> where u32: From<N> {
    fn n_nodes(&self) -> usize {
        self.row_ptrs.len().saturating_sub(1)
    }

    fn n_edges(&self) -> usize {
        self.column_indices.len()
    }

    /// Range of edge positions leaving `node`.
    fn row(&self, node: usize) -> (usize, usize) {
        (self.row_ptrs[node] as usize, self.row_ptrs[node + 1] as usize)
    }

    fn target(&self, edge: usize) -> usize {
        self.column_indices[edge] as usize
    }

    fn max_steps(&self, edge: usize) -> usize {
        (self.suffix_lens[edge] as usize).saturating_add(1)
    }

    /// Walks the path of `color` from the source to the sink.
    ///
    /// Returns the sequence length stored on the source edge and the NLZ
    /// dictionary as (colex rank, suffix length) pairs, in the order produced by
    /// [crate::encode_sequence].
    ///
    /// The returned value can be passed to [decode_path](super::decode_path) to recover the nucleotide
    /// sequence.
    ///
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError> {
        let mut sequence_length = 0_u32;
        let mut dictionary: Vec<(u32, u32)> = Vec::new();

        self.walk(choices, color, |node, pos| {
            if node == SOURCE_NODE as usize {
                sequence_length = self.suffix_lens[pos];
            } else {
                let colex_rank: u32 = self.node_weights[node].into();
                dictionary.push((colex_rank, self.suffix_lens[pos]));
            }
        })?;

        Ok((sequence_length, dictionary))
    }
}

