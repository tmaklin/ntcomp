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
use super::ColexGraphEdge;
use super::decode_sequence;

use std::cmp::Ordering;
use std::collections::{
    HashSet,
};

use blake3::Hash;

use petgraph::{
    Direction,
};
use petgraph::graph::{
    EdgeReference,
    Graph,
    NodeIndex,
};
use petgraph::visit::EdgeRef;

use sbwt::sbwt_index_variant::SbwtIndexVariant;

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
        visit_counts[0] = 1;
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
        self.node = edge.source();

        let visit_count = self.visit_counts[edge.target().index()];
        self.visited -= (visit_count == 1) as usize;
        self.visit_counts[edge.target().index()] = self.visit_counts[edge.target().index()].saturating_sub(1);
    }
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
