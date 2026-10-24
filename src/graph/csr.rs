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
use super::choices::Choices;
use super::paths::take_edge;
use super::{stored_entry, GraphDecodeError, GraphWalk, SINK_NODE, SOURCE_NODE};
use crate::nlz::{Dictionary, Entry};

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
    /// Returns the contig's dictionary (the sequence length stored on the
    /// source edge, and the entries in the order produced by
    /// [crate::nlz::parse]); [Dictionary::decode] recovers the sequence.
    ///
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<Dictionary, GraphDecodeError> {
        let mut sequence_len = 0_usize;
        let mut entries: Vec<Entry> = Vec::new();

        self.walk(choices, color, |node, pos| {
            if node == SOURCE_NODE as usize {
                sequence_len = self.suffix_lens[pos] as usize;
            } else {
                entries.push(stored_entry(self.node_weights[node].into(), self.suffix_lens[pos]));
            }
        })?;

        Ok(Dictionary { sequence_len, entries })
    }
}

