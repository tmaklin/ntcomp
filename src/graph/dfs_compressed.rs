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
    HashMap,
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

#[derive(Debug)]
struct ChainEdge {
    chain_weight: u32,
    visit_count: usize,
    nodes_in_chain: usize,
    chain_path: Vec<NodeIndex>,
}

fn compress_linear_chains(
    g: Graph<u32, ColexGraphEdge>,
    color: u32,
) -> Graph<u32, ChainEdge> {
    let mut compressed = petgraph::Graph::new();
    let mut node_mapping: HashMap<NodeIndex, NodeIndex> = HashMap::new();

    let source_node = NodeIndex::from(0);
    let sink_node = NodeIndex::from(1);
    node_mapping.insert(source_node, source_node);
    node_mapping.insert(sink_node, sink_node);

    let mut junctions = HashSet::new();
    junctions.insert(source_node);
    junctions.insert(sink_node);

    compressed.add_node(0_u32);
    compressed.add_node(1_u32);

    for node in g.node_indices() {
        let in_deg = g.edges_directed(node, Direction::Incoming).collect::<Vec<_>>().len();
        let out_deg = g.edges_directed(node, Direction::Outgoing).collect::<Vec<_>>().len();
        if in_deg != 1 || out_deg != 1 {
            junctions.insert(node);
        }
    }

    let get_or_create = |orig: NodeIndex, comp: &mut Graph<u32, ChainEdge>, map: &mut HashMap<NodeIndex, NodeIndex>| -> NodeIndex {
        *map.entry(orig).or_insert_with(|| comp.add_node(g[orig]))
    };

    for &junc in &junctions {
        let outgoing_edges: Vec<_> = g.edges_directed(junc, Direction::Outgoing).collect();

        for edge_ref in outgoing_edges {
            let mut chain_weight = edge_ref.weight().weight;
            let mut visit_count = edge_ref.weight().visit_counts[color as usize] as usize;
            let mut nodes_in_chain = 1;
            let mut chain_path = vec![edge_ref.target()];
            let mut curr = edge_ref.target();

            while !junctions.contains(&curr) {
                let next_edge = match g.edges_directed(curr, Direction::Outgoing).next() {
                    Some(e) => e,
                    None => break,
                };
                chain_weight += next_edge.weight().weight;
                visit_count += next_edge.weight().visit_counts[color as usize] as usize;
                nodes_in_chain += 1;
                chain_path.push(next_edge.target());
                curr = next_edge.target();
            }

            let new_source = get_or_create(junc, &mut compressed, &mut node_mapping);
            let new_target = get_or_create(curr, &mut compressed, &mut node_mapping);

            let chain_edge = ChainEdge {
                chain_weight,
                visit_count,
                nodes_in_chain,
                chain_path,
            };
            compressed.add_edge(new_source, new_target, chain_edge);
        }
    }

    compressed
}

#[derive(Clone, Debug)]
struct StackState {
    pub weight: u32,
    pub visited: usize,
    pub visit_counts: Vec<usize>,
    pub path: Vec<NodeIndex>,
    pub hash: Hash,
    pub node: NodeIndex,
    pub path_len: usize,
}

impl PartialEq for StackState {
    fn eq(&self, other: &Self) -> bool {
        let mut allowed = true;
        self.visit_counts.iter().zip(other.visit_counts.iter()).for_each(|(x, y)| {
            allowed &= x >= y // TODO the visits are undercounted for some reason
        });
            self.weight == other.weight &&
            self.visited == other.visited &&
            self.path.len() == other.path.len() &&
            self.hash == other.hash &&
            allowed &&
            // self.visit_counts == other.visit_counts &&
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
        graph: &Graph<u32, ChainEdge>,
        first_node: NodeIndex,
    ) -> Self {
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
            path_len: 1,
        }
    }

    pub fn for_path(
        graph: &Graph<u32, ChainEdge>,
        color: u32,
        hash: Hash,
        path_len: usize,
        first_node: NodeIndex,
        target_node: NodeIndex,
        target_weight: u32,
    ) -> Self {
        let mut visit_counts = vec![0_usize; graph.node_count()];
        let allowed_nodes: Vec<_> = graph
            .edge_references()
            .flat_map(|e| {
                visit_counts[e.target().index()] += e.weight().visit_count;
                [e.source(), e.target()]
            }).collect();
        visit_counts[first_node.index()] += 1;
        let nodes_in_path = HashSet::<NodeIndex>::from_iter(allowed_nodes).len();

        // Target length is stored in the weight for the first edge
        let first_edge: Vec<u32> = graph
            .edges_directed(first_node, Direction::Outgoing)
            .map(|e| e.weight().chain_weight).collect();
        assert!(first_edge.len() == 1);

        // Target weight is multiplied by 2 to account for the first edge's weight
        let total_weight = target_weight;
        StackState {
            weight: total_weight,
            visited: nodes_in_path,
            visit_counts,
            path: vec![NodeIndex::new(0); path_len],
            hash,
            node: target_node,
            path_len,
        }
    }

    pub fn advance(
        &mut self,
        edge: EdgeReference<ChainEdge>,
    ) {
        let nodes_added = edge.weight().nodes_in_chain;
        self.weight += edge.weight().chain_weight;
        self.path_len += nodes_added;

        let node = edge.target();
        self.node = node;
        self.path.extend(edge.weight().chain_path.iter());

        let visit_count = self.visit_counts[node.index()];
        self.visited += (visit_count == 0) as usize;
        self.visit_counts[node.index()] += nodes_added;
    }

    pub fn regress(
        &mut self,
        edge: EdgeReference<ChainEdge>,
    ) {
        let nodes_removed = edge.weight().nodes_in_chain;
        self.weight -= edge.weight().chain_weight;
        self.path_len = self.path_len.saturating_sub(nodes_removed);

        for _ in 0..edge.weight().chain_path.len() {
            if self.path.len() > 1 {
                self.path.pop();
            }
        }
        self.node = edge.source();

        let visit_count = self.visit_counts[edge.target().index()];
        self.visited -= (visit_count == nodes_removed) as usize;
        self.visit_counts[edge.target().index()] = self.visit_counts[edge.target().index()].saturating_sub(nodes_removed);
    }
}

pub fn search(
    graph: &Graph<u32, ColexGraphEdge>,
    color: u32,
    target_hash: Hash,
    sbwt: &SbwtIndexVariant,
    path_len: u32,
) -> Option<Vec<NodeIndex>> {

    let filtered_edges = graph.edge_references()
        .filter(|e| e.weight().colors.contains(&color));

    let mut target_weight: u32 = 0;
    let mut filtered_graph = graph.clone();
    filtered_graph.clear_edges();
    for e in filtered_edges {
        if e.source() == NodeIndex::from(0) {
            target_weight = e.weight().weight;
            let mut new_edge = e.weight().clone();
            new_edge.weight = 0;
            filtered_graph.add_edge(e.source(), e.target(), new_edge);
        } else {
            filtered_graph.add_edge(e.source(), e.target(), e.weight().clone());
        }
    }

    let compressed_graph = compress_linear_chains(filtered_graph, color);

    let first_node = NodeIndex::from(0);
    let target_node = NodeIndex::from(1);

    let target = StackState::for_path(
        &compressed_graph,
        color,
        target_hash,
        path_len.try_into().unwrap(),
        first_node,
        target_node,
        target_weight,
    );

    let mut tmp_stack = target.clone();
    tmp_stack.path = Vec::new();
    eprintln!("Target: {:?}", tmp_stack);

    stacker::grow(32 * 1024 * 1024, || {
        backtracking_search(
            &compressed_graph,
            &graph,
            sbwt,
            &target,
            &mut StackState::new(&compressed_graph, first_node),
        )
    })
}

fn backtracking_search(
    graph: &Graph<u32, ChainEdge>,
    original_graph: &Graph<u32, ColexGraphEdge>,
    sbwt: &SbwtIndexVariant,
    target: &StackState,
    stack: &mut StackState,
) -> Option<Vec<NodeIndex>> {

    let state = target.partial_cmp(stack).unwrap();
    if state == Ordering::Less {
        return None
    }

    if state == Ordering::Equal {
        let nucleotides = decode_sequence(original_graph, stack.path.as_slice(), sbwt);
        stack.hash = blake3::hash(&nucleotides);
        let mut tmp_stack = stack.clone();
        tmp_stack.path = Vec::new();
        eprintln!("Hashed: {:?}", tmp_stack);

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
        if let Some(path) = backtracking_search(graph, original_graph, sbwt, target, stack) {
            return Some(path)
        }
        stack.regress(edge);
    }

    None
}
