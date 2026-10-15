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

use crate::graph::ColexGraphEdge;

use egui::CentralPanel;
use egui_graphs::Graph;
use petgraph::prelude::StableGraph;

use petgraph::Direction;
use petgraph::prelude::*;
use std::collections::HashMap;

/// Compresses linear chains of nodes (in_degree == 1 && out_degree == 1)
/// into a single edge containing a vector of combined edge weights.
fn compress_linear_chains(
    g: &petgraph::Graph<u32, ColexGraphEdge>,
) -> petgraph::Graph<u32, Vec<ColexGraphEdge>> {
    let mut compressed = petgraph::Graph::new();
    let mut node_mapping: HashMap<NodeIndex, NodeIndex> = HashMap::new();

    // 1. Identify all junction nodes:
    // A node is a junction if it's a branch point, merge point, graph source, or graph sink.
    let mut junctions = std::collections::HashSet::new();
    for node in g.node_indices() {
        let in_deg = g.edges_directed(node, Direction::Incoming).collect::<Vec<_>>().len();
        let out_deg = g.edges_directed(node, Direction::Outgoing).collect::<Vec<_>>().len();
        if in_deg != 1 || out_deg != 1 {
            junctions.insert(node);
        }
    }

    // Helper closure to map old graph node indices to the new compressed graph indices
    let get_or_create = |orig: NodeIndex, comp: &mut petgraph::Graph<u32, Vec<ColexGraphEdge>>, map: &mut HashMap<NodeIndex, NodeIndex>| -> NodeIndex {
        *map.entry(orig).or_insert_with(|| comp.add_node(g[orig].clone()))
    };

    // 2. Traverse outward from each junction node to find linear paths
    for &junc in &junctions {
        // Collect outgoing edges so we release the borrow on `g` during traversal
        let outgoing_edges: Vec<_> = g.edges_directed(junc, Direction::Outgoing).collect();

        for edge_ref in outgoing_edges {
            let mut chain_weights = vec![edge_ref.weight().clone()];
            let mut curr = edge_ref.target();

            // Walk down the linear chain until we hit another junction node
            while !junctions.contains(&curr) {
                let next_edge = match g.edges_directed(curr, Direction::Outgoing).next() {
                    Some(e) => e,
                    None => break, // Reached a dead end
                };
                chain_weights.push(next_edge.weight().clone());
                curr = next_edge.target();
            }

            // Map both ends to the compressed graph and add the compressed edge
            let new_source = get_or_create(junc, &mut compressed, &mut node_mapping);
            let new_target = get_or_create(curr, &mut compressed, &mut node_mapping);

            compressed.add_edge(new_source, new_target, chain_weights);
        }
    }

    compressed
}

pub struct BasicApp {
    g: Graph<u32, Vec<ColexGraphEdge>>,
}

impl BasicApp {
    pub fn new(
        graph: petgraph::Graph<u32, ColexGraphEdge>,
        color: u32,
    ) -> Self {

        let filtered_edges = graph.edge_references()
                                  .filter(|e| e.weight().colors.contains(&color));

        let mut filtered_graph = graph.clone();
        filtered_graph.clear_edges();
        for e in filtered_edges {
            filtered_graph.add_edge(e.source(), e.target(), e.weight().clone());
        }


        let compressed = compress_linear_chains(&filtered_graph);
        let stable_g: StableGraph<u32, Vec<ColexGraphEdge>> = StableGraph::from(compressed);
        let g = Graph::from(&stable_g);
        Self { g }
    }
}

use egui_graphs::{GraphView, LayoutHierarchical, LayoutStateHierarchical};
type HLayout = LayoutHierarchical;
type HState = LayoutStateHierarchical;

impl eframe::App for BasicApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        CentralPanel::default().show(ui, |ui| {
            // Render the interactive graph view
            GraphView::<HState, HLayout>::new()
                .with_navigations(&egui_graphs::SettingsNavigation::default()
                                  .with_zoom_and_pan_enabled(true)
                                  .with_fit_to_screen_enabled(false))
                .with_interactions(&egui_graphs::SettingsInteraction::default()
                                   .with_dragging_enabled(true)
                                   .with_hover_enabled(true)
                                   .with_node_clicking_enabled(true)
                                   .with_node_selection_enabled(true)
                                   .with_node_selection_multi_enabled(true)
                                   .with_edge_clicking_enabled(true)
                                   .with_edge_selection_enabled(true)
                                   .with_edge_selection_multi_enabled(true)
                )
                .show(ui, &mut self.g);
        });
    }
}
