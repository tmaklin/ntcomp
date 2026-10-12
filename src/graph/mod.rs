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
//! NLZ-OG: NLZ dictionaries stored as coloured paths in a shared graph.
//!
//! ## Graph structure
//!
//! Node weights are colexicographic ranks of _k_-mers in the SBWT; node 0 is a
//! dummy source and node 1 a dummy sink. An edge `(u, v, s)` means "the
//! dictionary entry for `u` has suffix length `s` and is followed by `v`".
//! Edges are deduplicated across colours (contigs) and carry the set of colours
//! that use them.
//!
//! ## Decoding the graph
//!
//! The edge set alone does not determine the order in which a colour walks
//! the graph, so the encoder also stores a *choice stream*: whenever a colour
//! leaves a node that has more than one outgoing edge carrying that colour,
//! the rank of the taken edge among those edges is recorded. Decoding a
//! colour is a single linear walk from the source to the sink, and every
//! colour has its own stream so contigs can be decoded independently.
//!
use core::ops::Range;

use std::io::Write;

use indexmap::{IndexMap, IndexSet};

use petgraph::EdgeType;
use petgraph::graph::{
    Graph,
    NodeIndex,
};
use petgraph::visit::EdgeRef;

use sbwt::sbwt_index_variant::SbwtIndexVariant;

type E = Box<dyn std::error::Error>;

/// Node index of the dummy source node.
const SOURCE_NODE: u32 = 0;
/// Node index of the dummy sink node.
const SINK_NODE: u32 = 1;
/// Colex rank stored in the dummy source node. Real colex ranks must be smaller.
const SOURCE_COLEX: u32 = u32::MAX - 1;
/// Colex rank stored in the dummy sink node. Real colex ranks must be smaller.
const SINK_COLEX: u32 = u32::MAX;

#[derive(Clone, Debug, PartialEq, Eq, Hash, serde::Deserialize, serde::Serialize)]
pub struct ColexGraphEdge {
    pub weight: u32,
    pub colors: Vec<u32>,
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

/// Builds the coloured graph and remembers the path each colour takes.
pub struct GraphBuilder {

    /// colex rank -> node index (position in the set)
    nodes: IndexSet<u32>,

    /// (source node, target node, suffix length) -> sorted colours using the edge
    edges: IndexMap<(u32, u32, u32), Vec<u32>>,

    /// Per colour: edge ids (indexes into `edges`) in traversal order
    paths: Vec<Vec<usize>>,
}

impl Default for GraphBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl GraphBuilder {
    pub fn new() -> Self {
        let mut nodes = IndexSet::new();
        nodes.insert(SOURCE_COLEX);
        nodes.insert(SINK_COLEX);

        GraphBuilder {
            nodes,
            edges: IndexMap::new(),
            paths: Vec::new()
        }
    }

    fn insert_node(
        &mut self,
        colex_rank: usize
    ) -> u32 {
        let colex_rank: u32 = colex_rank.try_into().expect("colex rank does not fit in u32");
        assert!(colex_rank < SOURCE_COLEX, "colex rank collides with dummy node ranks");
        self.nodes.insert_full(colex_rank).0 as u32
    }

    fn edge(
        &mut self,
        source: u32,
        target: u32,
        weight: u32,
        color: u32,
    ) -> usize {
        let entry = self.edges.entry((source, target, weight));
        let id = entry.index();
        let colors = entry.or_default();
        // Colours are added in increasing order, so the colour sets stay sorted.
        if colors.last() != Some(&color) {
            colors.push(color);
        }
        id
    }

    /// Adds the dictionary of the next colour.
    ///
    /// `dictionary` is the output of [crate::encode_sequence]. The edge out of
    /// the source node stores `sequence_length`.
    pub fn add_path(
        &mut self,
        sequence_length: usize,
        dictionary: &[(usize, Range<usize>)],
    ) -> u32 {
        let color: u32 = self.paths.len().try_into().unwrap();
        let mut path = Vec::with_capacity(dictionary.len() + 1);

        let mut prev = SOURCE_NODE;
        let mut prev_weight: u32 = sequence_length.try_into().expect("sequence length does not fit in u32");
        for (suffix_len, colex) in dictionary {
            let node = self.insert_node(colex.start);
            path.push(self.edge(prev, node, prev_weight, color));
            prev = node;
            prev_weight = (*suffix_len).try_into().unwrap();
        }
        path.push(self.edge(prev, SINK_NODE, prev_weight, color));

        self.paths.push(path);
        color
    }

    /// Converts into the stored representation (graph + choice stream).
    pub fn to_csr(self) -> (Csr<u32>, Choices) {
        let n_nodes = self.nodes.len();
        let n_edges = self.edges.len();

        // Canonical edge order: by source, then target, then suffix length.
        let mut order: Vec<usize> = (0..n_edges).collect();
        order.sort_unstable_by_key(|&id| *self.edges.get_index(id).unwrap().0);

        let mut csr_pos = vec![0_u32; n_edges];
        let mut row_ptrs = vec![0_u32; n_nodes + 1];
        let mut column_indices = Vec::with_capacity(n_edges);
        let mut suffix_lens = Vec::with_capacity(n_edges);
        let mut colorsets = Vec::with_capacity(n_edges);

        for (pos, &id) in order.iter().enumerate() {
            let (&(source, target, weight), colors) = self.edges.get_index(id).unwrap();
            csr_pos[id] = pos as u32;
            row_ptrs[source as usize + 1] += 1;
            column_indices.push(target);
            suffix_lens.push(weight);
            colorsets.push(colors.clone());
        }
        for i in 0..n_nodes {
            row_ptrs[i + 1] += row_ptrs[i];
        }

        let csr = Csr {
            node_weights: self.nodes.into_iter().collect(),
            row_ptrs,
            column_indices,
            suffix_lens,
            colorsets,
        };

        // Record which edge each colour takes at every node where it has a choice.
        let mut streams: Vec<Vec<u32>> = Vec::with_capacity(self.paths.len());
        for (color, path) in self.paths.iter().enumerate() {
            let color = color as u32;
            let mut stream = Vec::new();
            for &id in path {
                let source = self.edges.get_index(id).unwrap().0.0 as usize;
                let pos = csr_pos[id] as usize;
                let mut n_options = 0;
                let mut rank = 0;
                for option in csr.outgoing_edges(source, color) {
                    if option < pos {
                        rank += 1;
                    }
                    n_options += 1;
                }
                if n_options > 1 {
                    stream.push(rank);
                }
            }
            streams.push(stream);
        }

        (csr, Choices { streams })
    }
}

/// Graph in compressed sparse row format; this is what gets stored.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct Csr<N> {
    pub node_weights: Vec<N>,
    pub row_ptrs: Vec<u32>,
    pub column_indices: Vec<u32>,
    pub suffix_lens: Vec<u32>,
    pub colorsets: Vec<Vec<u32>>,
}

impl<N> Csr<N> {
    /// Range of edge positions leaving `node`.
    pub fn row(
        &self,
        node: usize,
    ) -> (usize, usize) {
        (self.row_ptrs[node] as usize, self.row_ptrs[node + 1] as usize)
    }

    /// Edge positions leaving `node` that carry `color`.
    pub fn outgoing_edges(
        &self,
        node: usize,
        color: u32
    ) -> impl Iterator<Item = usize> + '_ {
        let (start, end) = self.row(node);
        (start..end).filter(move |&pos| self.colorsets[pos].binary_search(&color).is_ok())
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
            let (start, end) = self.row(source_idx);
            for pos in start..end {
                graph.add_edge(
                    NodeIndex::new(source_idx),
                    NodeIndex::new(self.column_indices[pos] as usize),
                    ColexGraphEdge {
                        weight: self.suffix_lens[pos],
                        colors: self.colorsets[pos].clone(),
                    },
                );
            }
        }

        graph
    }
}

/// Per-colour choice streams.
///
/// `streams[c][i]` is the choice made by colour `c` the _i_-th time it left a
/// node with more than one outgoing edge carrying `c`: the rank of the taken
/// edge among those edges, in CSR order.
///
/// Each colour has its own stream so any contig can be decoded on its own.
///
#[derive(serde::Serialize, serde::Deserialize, Debug, Default)]
pub struct Choices {
    pub streams: Vec<Vec<u32>>,
}

impl Choices {
    /// Total number of recorded choices.
    pub fn n_choices(&self) -> usize {
        self.streams.iter().map(|x| x.len()).sum()
    }
}

#[derive(Debug)]
pub struct GraphDecodeError(pub String);

impl std::fmt::Display for GraphDecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "graph decoding failed: {}", self.0)
    }
}

impl std::error::Error for GraphDecodeError {}

/// Walks the path of `color` from the source to the sink.
///
/// Returns the sequence length stored on the source edge and the NLZ
/// dictionary as (colex rank, suffix length) pairs, in the order produced by
/// [crate::encode_sequence].
///
/// The returned value can be passed to [decode_path] to recover the nucleotide
/// sequence.
///
/// `max_steps` bounds the walk (the stored path length) so a corrupt input
/// cannot loop forever.
///
pub fn walk_path(
    csr: &Csr<u32>,
    choices: &Choices,
    color: u32,
    max_steps: usize,
) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError> {
    let mut stream = choices
        .streams
        .get(color as usize)
        .ok_or_else(|| GraphDecodeError(format!("no choice stream for colour {}", color)))?
        .iter();
    let mut node = SOURCE_NODE as usize;
    let mut sequence_length = 0_u32;
    let mut dictionary: Vec<(u32, u32)> = Vec::with_capacity(max_steps);

    for _ in 0..max_steps {
        let mut options = csr.outgoing_edges(node, color);
        let first = options.next().ok_or_else(|| GraphDecodeError(format!("colour {} has no edge out of node {}", color, node)))?;
        let pos = match options.next() {
            None => first,
            Some(second) => {
                let rank = *stream.next().ok_or_else(|| GraphDecodeError(format!("choice stream of colour {} exhausted at node {}", color, node)))? as usize;
                match rank {
                    0 => first,
                    1 => second,
                    _ => options.nth(rank - 2).ok_or_else(|| GraphDecodeError(format!("invalid choice {} for colour {} at node {}", rank, color, node)))?,
                }
            },
        };

        if node == SOURCE_NODE as usize {
            sequence_length = csr.suffix_lens[pos];
        } else {
            dictionary.push((csr.node_weights[node], csr.suffix_lens[pos]));
        }

        node = csr.column_indices[pos] as usize;
        if node == SINK_NODE as usize {
            if stream.next().is_some() {
                return Err(GraphDecodeError(format!("unused choices left for colour {}", color)))
            }
            return Ok((sequence_length, dictionary))
        }
    }

    Err(GraphDecodeError(format!("colour {} did not reach the sink in {} steps", color, max_steps)))
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
