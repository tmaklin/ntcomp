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

use std::collections::{
    HashMap,
    HashSet,
};
use std::io::Write;

use indexmap::IndexSet;

use petgraph::graph::{
    Graph,
    NodeIndex,
};
use petgraph::visit::EdgeRef;

use sbwt::sbwt_index_variant::SbwtIndexVariant;

type E = Box<dyn std::error::Error>;

#[derive(Clone, Copy, Debug, serde::Serialize)]
pub struct ColexGraphEdge {
    pub weight: u32,
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

pub fn decode_path(
    path: Vec<(u32, u32)>,
    sbwt: &SbwtIndexVariant,
) -> Vec<u8> {
    let mut sequence: Vec<u8> = Vec::new();
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            path.into_iter().rev().for_each(|(colex_rank, suffix_len)| {
                let kmer = if suffix_len > k as u32 {
                    let kmer = sbwt.access_kmer(colex_rank as usize);
                    let new_kmer = crate::left_extend_kmer2(&kmer, sbwt, (suffix_len - k as u32) as usize);
                    assert_eq!(new_kmer.len(), suffix_len as usize);
                    new_kmer
                } else {
                    sbwt.access_kmer(colex_rank as usize)
                };
                sequence.extend(kmer[(kmer.len() - (suffix_len as usize))..kmer.len()].iter());
            });
        },
    };

    sequence
}

pub fn write_to<W: Write>(
    graph: &Graph<u32, ColexGraphEdge>,
    writer: &mut W,
) -> Result<(), E> {
    let bytes = postcard::to_allocvec(&graph)?;
    writer.write_all(&bytes)?;
    Ok(())
}

pub fn encode_to<W: Write>(
    graph: Graph<u32, ColexGraphEdge>,
    writer: &mut W,
) -> Result<(), E> {

    let mut deltas: Vec<i32> = Vec::with_capacity(graph.edge_count());
    let mut weights: Vec<u32> = Vec::with_capacity(graph.edge_count());
    let mut color_counts: HashMap<u32, u32> = HashMap::new();

    for e in graph.edge_references() {
        let from_index: i64 = e.source().index().try_into()?;
        let to_index: i64 = e.target().index().try_into()?;
        let delta: i32 = (from_index - to_index).try_into()?;
        deltas.push(delta);

        let color: u32 = e.weight().color;
        if color_counts.contains_key(&color) {
            *color_counts.get_mut(&color).unwrap() += 1;
        } else {
            color_counts.insert(color, 1_u32);
        }

        weights.push(e.weight().weight);
    }

    eprintln!("Max weight: {}", weights.iter().max().unwrap());
    // Nodes: encode the count as the node index is always just incremented by 1
    let node_count = graph.node_count();
    writer.write_all(&postcard::to_allocvec(&node_count)?)?;

    // Sources and targets: encode as difference?
    // TODO This requires storing the start index for every contig so we can start using the diffs
    writer.write_all(&postcard::to_allocvec(&deltas)?)?;

    // Colors: encode as (count, colour) since these are always contiguous
    let color_rle: Vec<(u32, u32)> = color_counts.into_iter().collect();
    writer.write_all(&postcard::to_allocvec(&color_rle)?)?;

    // Weights: use lzma
    writer.write_all(&postcard::to_allocvec(&weights)?)?;

    Ok(())
}
