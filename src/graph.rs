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
use core::ops::Range;

use std::collections::HashSet;
use std::io::Write;

use indexmap::IndexSet;

use petgraph::graph::{
    Graph,
    NodeIndex,
};
use petgraph::visit::EdgeRef;

type E = Box<dyn std::error::Error>;

#[derive(Debug, serde::Serialize)]
pub struct ColexGraphEdge {
    weight: u32,
    color: u32,
}

pub fn extract_path(
    graph: &Graph<u32, ColexGraphEdge>,
    want_color: u32,
) -> Vec<(u32, u32)> {
    graph.edge_references().filter(|x| x.weight().color == want_color).map(|edge| {
        let remapped_colex = graph[edge.source()];
        let match_length = edge.weight().weight;
        (remapped_colex, match_length)
    }).collect::<Vec<(u32, u32)>>()
}

pub fn insert_edge(
    graph: &mut Graph<u32, ColexGraphEdge>,
    node_indexes: &mut HashSet<u32>,
    colex_remapping: &mut IndexSet<u32>,
    color: u32,
    entry_from: &(usize, Range<usize>),
    entry_to: &(usize, Range<usize>),
) {
    let prev_idx: u32 = entry_from.1.start.try_into().unwrap();
    let curr_idx: u32 = entry_to.1.start.try_into().unwrap();
    let weight = entry_from.0;

    colex_remapping.insert(prev_idx);
    colex_remapping.insert(curr_idx);

    let edge_start: u32 = colex_remapping.get_index_of(&prev_idx).unwrap().try_into().unwrap();
    let edge_end: u32 = colex_remapping.get_index_of(&curr_idx).unwrap().try_into().unwrap();

    let from = if node_indexes.contains(&edge_start) {
        NodeIndex::from(edge_start)
    } else {
        node_indexes.insert(edge_start);
        graph.add_node(edge_start)
    };

    let to = if node_indexes.contains(&edge_end) {
        NodeIndex::from(edge_end)
    } else {
        node_indexes.insert(edge_end);
        graph.add_node(edge_end)
    };

    graph.try_add_edge(from, to, ColexGraphEdge { weight: weight.try_into().unwrap(), color }).unwrap();
}

pub fn write_to<W: Write>(
    graph: &Graph<u32, ColexGraphEdge>,
    writer: &mut W,
) -> Result<(), E> {
    let bytes = postcard::to_allocvec(&graph)?;
    writer.write_all(&bytes)?;
    Ok(())
}
