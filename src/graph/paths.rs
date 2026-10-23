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
//! Walking a colour's path through a stored graph with its choice stream.
//!
use super::choices::{ChoiceReader, Choices};
use super::GraphDecodeError;

/// A graph of nodes with ranked outgoing edges whose colours' paths can be
/// walked from the source to the sink: random access to any colour (contig)
/// stored in it. For the chain graph the nodes are the junctions and the
/// edges the chains.
///
/// Every dictionary entry covers at least one base, so the path of a contig of
/// `n` bases has at most `n + 1` edges (counting the edge out of the source,
/// which stores `n`). Walks stop with an error after that many edges, so a
/// corrupt input cannot loop forever.
pub trait GraphWalk {
    /// Number of nodes.
    fn n_nodes(&self) -> usize;

    /// Number of edges.
    fn n_edges(&self) -> usize;

    /// Edges `start..end` leaving `node`, in rank order.
    fn row(
        &self,
        node: usize,
    ) -> (usize, usize);

    /// The node `edge` leads to.
    fn target(
        &self,
        edge: usize,
    ) -> usize;

    /// Most steps a path that starts with `edge` can take: its contig length
    /// (the suffix length stored on the edges out of the source) + 1.
    fn max_steps(
        &self,
        edge: usize,
    ) -> usize;

    /// Walk the path of `color` from the source to the sink.
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError>;
}

/// Picks the edge out of `node`, whose outgoing edges are `start..end`,
/// reading the choice stream if there is more than one.
pub(super) fn take_edge(
    (start, end): (usize, usize),
    choices: &mut ChoiceReader<'_>,
    color: u32,
    node: usize,
) -> Result<usize, GraphDecodeError> {
    match end - start {
        0 => Err(GraphDecodeError(format!("colour {} has no edge out of node {}", color, node))),
        1 => Ok(start),
        n_options => {
            let rank = choices.next(n_options as u32).ok_or_else(|| GraphDecodeError(format!("choice stream of colour {} exhausted at node {}", color, node)))?;
            Ok(start + rank as usize)
        },
    }
}

