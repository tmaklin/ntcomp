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

use std::cmp::Ordering;
use std::collections::{
    HashMap,
    HashSet,
};
use std::io::Write;

use blake3::Hash;

use indexmap::IndexSet;

use petgraph::{
    Direction,
    EdgeType,
};
use petgraph::graph::{
    EdgeReference,
    Graph,
    NodeIndex,
};
use petgraph::visit::EdgeRef;

use sbwt::sbwt_index_variant::SbwtIndexVariant;

type E = Box<dyn std::error::Error>;

#[derive(Clone, Debug)]
struct StackState {
    pub weight: u32,
    pub visited: usize,
    pub visit_counts: Vec<usize>,
    pub path: Vec<NodeIndex>,
    pub hash: Hash,
    pub node: NodeIndex,
}

impl PartialEq for StackState {
    fn eq(&self, other: &Self) -> bool {
            self.weight == other.weight &&
            self.visited == other.visited &&
            self.path.len() == other.path.len() &&
            self.hash == other.hash &&
            self.visit_counts == other.visit_counts &&
            self.node == other.node
    }
}


impl PartialOrd for StackState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let is_less =
            self.weight < other.weight ||
            self.visited < other.visited ||
            self.path.len() < other.path.len();

        let is_greater =
            self.weight > other.weight ||
            self.visited > other.visited ||
            self.path.len() > other.path.len();

        let is_partial_eq =
            self.weight == other.weight &&
            self.visited == other.visited &&
            self.path.len() == other.path.len() &&
            self.node == other.node;

        if is_partial_eq {
            Some(Ordering::Equal)
        } else if is_less {
            Some(Ordering::Less)
        } else if is_greater {
            Some(Ordering::Greater)
        } else {
            None
        }
    }
}

impl StackState {
    pub fn new(
        graph: &Graph<u32, ColexGraphEdge>,
    ) -> Self {
        let first_node = NodeIndex::from(0_u32);
        let mut path = Vec::with_capacity(graph.node_count());
        path.push(first_node);
        let mut visit_counts = vec![0_usize; graph.node_count()];
        visit_counts[first_node.index()] = 1;

        StackState {
            weight: 0,
            visited: 1,
            visit_counts,
            path,
            hash: Hash::from_bytes([0_u8; 32]),
            node: first_node,
        }
    }

    pub fn for_path(
        graph: &Graph<u32, ColexGraphEdge>,
        color: u32,
        hash: Hash,
        path_len: usize,
    ) -> Self {
        let mut visit_counts = vec![0_usize; graph.node_count()];
        let allowed_edges: Vec<_> = graph
            .edge_references()
            .filter(|e| e.weight().colors.contains(&color))
            .flat_map(|e| {
                visit_counts[e.target().index()]  = e.weight().visit_counts[color as usize] as usize;
                [e.source(), e.target()]
            }).collect();
        visit_counts[0] += 1;
        let nodes_in_path = HashSet::<NodeIndex>::from_iter(allowed_edges).len();

        // Paths for all colors start at dummy NodeIndex 0
        let first_node = NodeIndex::from(0_u32);

        // Target length is stored in the weight for the first edge
        let first_edge: Vec<u32> = graph.edges_directed(first_node, Direction::Outgoing)
                                        .filter(|e| e.weight().colors.contains(&color))
                                        .map(|e| e.weight().weight)
                                        .collect();
        assert!(first_edge.len() == 1);

        // Target weight is multiplied by 2 to account for the first edge's weight
        let total_weight = first_edge[0] * 2;

        // Paths for all colors end at dummy NodeIndex 1
        let target_node = NodeIndex::from(1_u32);
        StackState {
            weight: total_weight,
            visited: nodes_in_path,
            visit_counts,
            path: vec![NodeIndex::new(0); path_len],
            hash,
            node: target_node,
        }
    }

    pub fn advance(
        &mut self,
        edge: EdgeReference<ColexGraphEdge>,
    ) {
        self.weight += edge.weight().weight;

        let node = edge.target();
        self.node = node;
        self.path.push(node);

        let visit_count = self.visit_counts[node.index()];
        self.visited += (visit_count == 0) as usize;
        self.visit_counts[node.index()] += 1;
    }

    pub fn regress(
        &mut self,
        edge: EdgeReference<ColexGraphEdge>,
    ) {
        self.weight -= edge.weight().weight;

        if self.path.len() > 1 {
            self.path.pop();
        }
        let prev = edge.source();
        self.node = prev;

        let visit_count = self.visit_counts[prev.index()];
        self.visited -= (visit_count == 0) as usize;
        self.visit_counts[prev.index()] = self.visit_counts[prev.index()].saturating_sub(1);
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, serde::Deserialize, serde::Serialize)]
pub struct ColexGraphEdge {
    weight: u32,
    colors: Vec<u32>,
    visit_counts: Vec<u32>,
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
    for i in 1..nodes.len() {
        let node_idx = nodes[i];
        let colex_rank = graph[node_idx];
        let edges = graph.edges_directed(node_idx, Direction::Outgoing);
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
    target_hash: Hash,
    sbwt: &SbwtIndexVariant,
    path_len: u32,
) -> Option<Vec<NodeIndex>> {

    let target = StackState::for_path(
        graph,
        color,
        target_hash,
        path_len.try_into().unwrap(),
    );

    let filtered_edges = graph.edge_references()
        .filter(|e| e.weight().colors.contains(&color));

    let mut filtered_graph = graph.clone();
    filtered_graph.clear_edges();
    for e in filtered_edges {
        filtered_graph.add_edge(e.source(), e.target(), e.weight().clone());
    }

    stacker::grow(32 * 1024 * 1024, || {
        backtracking_search(
            &filtered_graph,
            sbwt,
            &target,
            &mut StackState::new(graph),
        )
    })
}

fn backtracking_search(
    graph: &Graph<u32, ColexGraphEdge>,
    sbwt: &SbwtIndexVariant,
    target: &StackState,
    stack: &mut StackState,
) -> Option<Vec<NodeIndex>> {
    let state = target.partial_cmp(stack).unwrap();
    if state == Ordering::Less {
        return None
    }

    if state == Ordering::Equal {
        let nucleotides = decode_sequence(graph, stack.path.as_slice(), sbwt);
        stack.hash = blake3::hash(&nucleotides);
        if target.eq(stack) {
            return Some(std::mem::take(&mut stack.path))
        } else {
            return None
        }
    }

    let outgoing_edges: Vec<_> = graph
        .edges_directed(stack.node, petgraph::Direction::Outgoing)
        .filter(|e| {
            let target_idx = e.target().index();
            let visit_count: usize = stack.visit_counts[target_idx];
            let allowed_count: usize = target.visit_counts[target_idx];
            visit_count < allowed_count
        }).collect();

    for edge in outgoing_edges {
        stack.advance(edge);
        if let Some(path) = backtracking_search(graph, sbwt, target, stack) {
            return Some(path)
        }
        stack.regress(edge);
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

    graph.try_add_edge(from, to, ColexGraphEdge { weight: suffix_len, colors: vec![color], visit_counts: vec![1] }).unwrap();
}

pub fn deduplicate_edges(
    graph: &mut Graph<u32, ColexGraphEdge>,
) -> u32 {

    let mut visit_counts: HashMap<NodeIndex, Vec<u32>> = HashMap::new();

    let first_node = NodeIndex::from(0_u32);
    let n_colors: usize = HashSet::<u32>::from_iter(graph.edges_directed(first_node, Direction::Outgoing).flat_map(|e| e.weight().colors.clone()).collect::<Vec<u32>>()).len();

    let mut edge_colors: HashMap<(NodeIndex, NodeIndex, u32), Vec<u32>> = HashMap::from_iter(
        graph.edge_references().map(|x| {
            x.weight().colors.iter().for_each(|color| {
                visit_counts.entry(x.target()).and_modify(|e| e[*color as usize] += 1).or_insert( {
                    let mut counts = vec![0; n_colors];
                    counts[*color as usize] = 1;
                    counts
                });
            });
            ((x.source(), x.target(), x.weight().weight), x.weight().colors.clone())
        })
    );

    let max_visits = *visit_counts.values().flatten().max().unwrap();

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
        let visit_counts = visit_counts.get(&x.target()).unwrap().clone();
        (x.source(),
         x.target(),
         ColexGraphEdge {
             weight: x.weight().weight,
             colors,
             visit_counts,
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

/// CSR for storage
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct Csr<N> {
    pub node_weights: Vec<N>,
    pub row_ptrs: Vec<u32>,
    pub column_indices: Vec<u32>,
    pub suffix_lens: Vec<u32>,
    pub colorsets: Vec<Vec<u32>>,
    pub visit_counts: Vec<Vec<u32>>,
}

impl<N> Csr<N> {
    pub fn from_petgraph<Ty: EdgeType>(
        graph: Graph<N, ColexGraphEdge, Ty>,
    ) -> Self
    where
        N: Copy,
    {
        let node_count = graph.node_count();
        let edge_count = graph.edge_count();

        let mut row_ptrs = Vec::with_capacity(node_count + 1);
        let mut column_indices = Vec::with_capacity(edge_count);
        let mut suffix_lens = Vec::with_capacity(edge_count);
        let mut node_weights = Vec::with_capacity(node_count);
        let mut colorsets: Vec<Vec<u32>> = Vec::with_capacity(edge_count);
        let mut visit_counts: Vec<Vec<u32>> = Vec::with_capacity(edge_count);

        let mut current_offset = 0;
        row_ptrs.push(current_offset);

        for node in graph.node_indices() {
            let edges = graph.edges_directed(node, Direction::Outgoing);

            for edge in edges {
                column_indices.push(edge.target().index().try_into().unwrap());
                suffix_lens.push(edge.weight().weight);
                colorsets.push(edge.weight().colors.clone());
                visit_counts.push(edge.weight().visit_counts.clone());
                current_offset += 1;
            }
            node_weights.push(graph[node]);
            row_ptrs.push(current_offset);
        }

        Self {
            node_weights,
            row_ptrs,
            column_indices,
            suffix_lens,
            colorsets,
            visit_counts,
        }
    }

    pub fn to_petgraph<Ty: EdgeType>(
        self,
    ) -> Graph<N, ColexGraphEdge, Ty>
    where
        N: Copy,
    {
        let mut graph = Graph::with_capacity(self.node_weights.len(), self.column_indices.len());

        for weight in &self.node_weights {
            graph.add_node(*weight);
        }

        for source_idx in 0..self.node_weights.len() {
            let start = self.row_ptrs[source_idx];
            let end = self.row_ptrs[source_idx + 1];

            for edge_offset in start..end {
                let target_idx = self.column_indices[edge_offset as usize] as usize;
                let suffix_len = &self.suffix_lens[edge_offset as usize];
                let colorset = &self.colorsets[edge_offset as usize];
                let visit_counts = &self.visit_counts[edge_offset as usize];

                graph.add_edge(
                    NodeIndex::new(source_idx),
                    NodeIndex::new(target_idx),
                    ColexGraphEdge {
                        weight: *suffix_len,
                        colors: colorset.to_vec(),
                        visit_counts: visit_counts.clone(),
                    },
                );
            }
        }

        graph
    }
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
