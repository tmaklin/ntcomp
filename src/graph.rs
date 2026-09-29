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

use blake3::Hash;

use indexmap::IndexSet;

use petgraph::graph::{
    Graph,
    NodeIndex,
};
use petgraph::visit::EdgeRef;

use sbwt::sbwt_index_variant::SbwtIndexVariant;

type E = Box<dyn std::error::Error>;

#[derive(Clone, Debug, PartialEq, Eq, Hash, serde::Deserialize, serde::Serialize)]
pub struct ColexGraphEdge {
    pub weight: u32,
    colors: Vec<u32>,
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

pub fn decode_sequence(
    graph: &Graph<u32, ColexGraphEdge>,
    nodes: &[NodeIndex],
    sbwt: &SbwtIndexVariant,
) -> Vec<u8> {
    let mut path: Vec<(u32, u32)> = Vec::new();
    for i in 0..nodes.len() {
        let node_idx = nodes[i];
        let colex_rank = graph[node_idx];
        let edges = graph.edges_directed(node_idx, petgraph::Direction::Outgoing);
        let tmp = if i < nodes.len() - 1 { i + 1 } else { 0 };
        let next_node_idx = nodes[tmp];
        for e in edges {
            if e.target() == next_node_idx {
                let suffix_len: u32 = e.weight().weight;
                path.push((colex_rank, suffix_len));
                break;
            }
        }
    }

    decode_path(path, sbwt)
}

pub fn search(
    graph: &Graph<u32, ColexGraphEdge>,
    color: u32,
    total_weight: u32,
    hash: Hash,
    sbwt: &SbwtIndexVariant,
    max_visits: u32,
) -> Option<Vec<NodeIndex>> {

    let nodes_in_path = HashSet::<NodeIndex>::from_iter(graph
        .edge_references()
        .filter(|e| e.weight().colors.contains(&color))
        .flat_map(|e| [e.source(), e.target()])
        .collect::<Vec<NodeIndex>>()
    ).len();

    // Dummy nodes that denote start/end
    let first_node = NodeIndex::from(0_u32);
    let last_node = NodeIndex::from(1_u32);

    stacker::grow(1024 * 1024 * 1024, || {
        backtracking_search(
            graph,
            first_node,
            last_node,
            color,
            total_weight,
            nodes_in_path,
            hash,
            &mut HashMap::new(),
            &mut 0_u32,
            &mut Vec::new(),
            sbwt,
            max_visits as usize,
        )
    })
}

pub fn backtracking_search(
    graph: &Graph<u32, ColexGraphEdge>,
    current: NodeIndex,
    end: NodeIndex,
    color: u32,
    target_weight: u32,
    nodes_to_visit: usize,
    hash: Hash,
    visit_counts: &mut HashMap<NodeIndex, usize>,
    current_weight: &mut u32,
    path: &mut Vec<NodeIndex>,
    sbwt: &SbwtIndexVariant,
    max_visits: usize,
) -> Option<Vec<NodeIndex>> {
    if *current_weight > target_weight || visit_counts.len() > nodes_to_visit {
        return None
    }

    path.push(current);
    visit_counts.entry(current).and_modify(|e| *e += 1).or_insert(1);

    if current == end && *current_weight == target_weight {
        let nucleotides = decode_sequence(graph, &path.clone(), sbwt);
        let hash_got = blake3::hash(&nucleotides);
        if hash_got == hash {
            return Some(path.to_vec())
        } else {
            return None
        }
    }

    let outgoing_edges = graph
        .edges_directed(current, petgraph::Direction::Outgoing)
        .filter(|e| e.weight().colors.contains(&color));

    for edge in outgoing_edges {
        *current_weight += edge.weight().weight;

        let is_valid = visit_counts.get(&edge.target()).unwrap_or(&0_usize) <= &max_visits;

        if is_valid {
            if backtracking_search(graph, edge.target(), end, color, target_weight, nodes_to_visit, hash, visit_counts, current_weight, path, sbwt, max_visits).is_some() {
                return Some(path.to_vec())
            }
        }

        *current_weight -= edge.weight().weight;
    }

    path.pop();
    let count = *visit_counts.get(&current).unwrap();
    if count == 1 {
        visit_counts.remove(&current);
    } else {
        *visit_counts.get_mut(&current).unwrap() -= 1;
    }

    None
}

pub fn insert_edge(
    graph: &mut Graph<u32, ColexGraphEdge>,
    node_indexes: &mut IndexSet<u32>,
    color: u32,
    entry_from: &(usize, Range<usize>),
    entry_to: &(usize, Range<usize>),
) {
    // TODO should check that we don't use colex ranks 0 or 1,
    // iirc these are used for dummy nodes in the SBWT so should be fine
    let source_colex: u32 = entry_from.1.start.try_into().unwrap();
    let target_colex: u32 = entry_to.1.start.try_into().unwrap();

    let suffix_len: u32 = entry_from.0.try_into().unwrap();

    let from: NodeIndex<u32> = if node_indexes.contains(&source_colex) {
        let source_index = node_indexes.get_index_of(&source_colex).unwrap();
        NodeIndex::from(source_index as u32)
    } else {
        node_indexes.insert(source_colex);
        graph.add_node(source_colex)
    };

    let to = if node_indexes.contains(&target_colex) {
        let target_index = node_indexes.get_index_of(&target_colex).unwrap();
        NodeIndex::from(target_index as u32)
    } else {
        node_indexes.insert(target_colex);
        graph.add_node(target_colex)
    };

    graph.try_add_edge(from, to, ColexGraphEdge { weight: suffix_len, colors: vec![color] }).unwrap();
}

pub fn deduplicate_edges(
    graph: &mut Graph<u32, ColexGraphEdge>,
) -> u32 {

    let mut visit_counts: HashMap<(NodeIndex, NodeIndex, u32), u32> = HashMap::new();

    let mut edge_colors: HashMap<(NodeIndex, NodeIndex, u32), Vec<u32>> = HashMap::from_iter(
        graph.edge_references().map(|x| {
            let key = (x.source(), x.target(), x.weight().weight);
            visit_counts.entry(key).and_modify(|e| *e += 1).or_insert(1);
            ((x.source(), x.target(), x.weight().weight), x.weight().colors.clone())
        })
    );

    let max_visits = visit_counts.into_iter().map(|(_, val)| val).max().unwrap();

    for e in graph.edge_references() {
        let key = (e.source(), e.target(), e.weight().weight);
        if !edge_colors.get(&key).unwrap().contains(&e.weight().colors[0]) {
            edge_colors.get_mut(&key).unwrap().push(e.weight().colors[0]);
            edge_colors.get_mut(&key).unwrap().sort();
        }
    }

    let edges_iter = graph.edge_references().map(|x| {
        let key = (x.source(), x.target(), x.weight().weight);
        let colors = edge_colors.get(&key).unwrap().clone();
        (x.source(),
         x.target(),
         ColexGraphEdge {
             weight: x.weight().weight,
             colors,
         },
        )
    });

    let unique_edges: IndexSet<(NodeIndex, NodeIndex, ColexGraphEdge)> = IndexSet::from_iter(edges_iter);

    graph.clear_edges();
    for (source, target, weight) in unique_edges {
        graph.try_add_edge(source, target, weight).unwrap();
    }

    max_visits
}

/// With colex ranks stored in the node
pub fn encode_to<W: Write>(
    graph: &Graph<u32, ColexGraphEdge>,
    colex_remapping: &IndexSet<u32>,
    writer: &mut W,
) -> Result<(), E> {

    let mut deltas: Vec<i32> = Vec::with_capacity(graph.edge_count());
    let mut weights: Vec<u32> = Vec::with_capacity(graph.edge_count());
    let mut colors: Vec<Vec<u32>> = Vec::with_capacity(graph.edge_count());

    for e in graph.edge_references() {
        let from_index: usize = e.source().index();
        let to_index: usize = e.target().index();
        let source_colex: i64 = (*colex_remapping.get_index(from_index).unwrap()).into();
        let target_colex: i64 = (*colex_remapping.get_index(to_index).unwrap()).into();
        let delta: i32 = (target_colex - source_colex).try_into()?;
        deltas.push(delta);

        colors.push(e.weight().colors.clone());

        let suffix_len = e.weight().weight;
        weights.push(suffix_len);
    }

    deltas.shrink_to_fit();
    weights.shrink_to_fit();
    colors.shrink_to_fit();
    // // Nodes: encode the count as the node index is always just incremented by 1
    // let node_count = graph.node_count();
    // writer.write_all(&postcard::to_allocvec(&node_count)?)?;

    // // Colors: encode as (count, colour) since these are always contiguous
    // let color_rle: Vec<(u32, u32)> = color_counts.into_iter().collect();
    // writer.write_all(&postcard::to_allocvec(&color_rle)?)?;

    // Sources and targets: encode as difference?
    // TODO This requires storing the start index for every contig so we can start using the diffs
    writer.write_all(&postcard::to_allocvec(&deltas)?)?;

    // Weights: use lzma
    writer.write_all(&postcard::to_allocvec(&weights)?)?;

    // Colors
    writer.write_all(&postcard::to_allocvec(&colors)?)?;

    // writer.write_all(&postcard::to_allocvec(&paths)?)?;
    Ok(())
}
