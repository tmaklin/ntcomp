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
//! Node weights are colexicographic ranks of _k_-mers in the SBWT. Node 0 is a
//! dummy source and node 1 a dummy sink. An edge `(u, v, s)` means "the
//! dictionary entry for `u` has suffix length `s` and is followed by `v`".
//! Edges are deduplicated across colours (contigs). Colours used by the edges
//! are not stored but can be retrieved by walking the paths.
//!
//! ## Decoding the graph
//!
//! The edge set alone does not determine the path a colour takes through the
//! graph, so the encoder also stores a *choice stream*: whenever a colour
//! leaves a node that has more than one outgoing edge, the rank of the taken
//! edge among all of them is recorded. The edges of every node are ordered by
//! how many times the colours take them, so the most common choice has rank 0
//! and the ranks can be compressed, see [choices].
//!
//! Decoding a colour is a single linear walk from the source to the sink, and
//! every colour has its own stream so contigs can be decoded independently. The
//! colours using each edge can be recovered by walking every colour, see
//! [Csr::colorsets].
//!
//! ## Collapsed chains
//!
//! Optionally, linear chains of nodes are stored as single edges between
//! junctions, see [collapsed]. Linear chains are paths on the graph composed
//! entirely of inner nodes, meaning nodes that have exactly one incoming and
//! exactly one outgoing edge. All other nodes are junction nodes.
//!
use core::ops::Range;

use std::io::Write;

use indexmap::IndexSet;

use petgraph::EdgeType;
use petgraph::graph::{
    Graph,
    NodeIndex,
};
use petgraph::visit::EdgeRef;

use sbwt::sbwt_index_variant::SbwtIndexVariant;

pub mod bits;
pub mod choices;
pub mod collapsed;
pub mod decoder;
pub mod encoder;
pub mod huffman;
pub mod lzma_wrapper;

pub use choices::Choices;
use choices::{ChoiceReader, ChoiceWriter};
use collapsed::ChainCsr;

type E = Box<dyn std::error::Error>;

/// Node index of the dummy source node.
const SOURCE_NODE: u32 = 0;
/// Node index of the dummy sink node.
const SINK_NODE: u32 = 1;
/// Colex rank stored in the dummy source node. Real colex ranks must be smaller.
const SOURCE_COLEX: u32 = u32::MAX - 1;
/// Colex rank stored in the dummy sink node. Real colex ranks must be smaller.
const SINK_COLEX: u32 = u32::MAX;

/// [Implicit](SuffixLen::Implicit) suffix lengths are stored as a 0; a real suffix length can never be 0.
const IMPLICIT_SUFFIX_LEN: u32 = 0;

/// Walkable graph that implements random access to any colour (contig) stored inside.
///
/// Every dictionary entry covers at least one base, so the path of a contig of
/// `n` bases has at most `n + 1` edges (counting the edge out of the source,
/// which stores `n`). Walks stop with an error after that many edges, so a
/// corrupt input cannot loop forever.
pub trait Walkable {
    /// Walk the path of `color` from the source to the sink.
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError>;
}

/// Suffix length and type for an edge
///
/// Classifies the suffix length for a dictionary entry by specifying whether
/// its implicit (for a left extended _k_-mer but can be inferred), explicit
/// (left extended and needs to be declared), or plain (subsequence of a plan
/// _k_-mer). Stores the value for explicit and plain types.
///
pub enum SuffixLen {
    /// The dictionary entry is a left extended _k_-mer whose extension stopped
    /// where [crate::left_extend_kmer2] stops, so [decode_path] recovers the
    /// length by extending the _k_-mer as far as possible.
    Implicit,

    /// Entry is a left extended _k_-mer, whose extension stops after an
    /// explicit number of extensions, so [decode_path] needs to know the
    /// length.
    Explicit { suffix_len: u32 },

    /// Entry is a plain _k_-mer whose suffix length is less or equal to the
    /// value of k.
    Plain { suffix_len: u32 },
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, serde::Deserialize, serde::Serialize)]
pub struct ColexGraphEdge {
    pub weight: u32,
    pub colors: Vec<u32>,
}

/// Decodes the nucleotide sequence of length `sequence_length` from its NLZ
/// dictionary (colex rank, suffix length), as returned by [walk_path].
///
/// The dictionary is processed in its own order, from the end of the sequence
/// to the start, so that the end position of every entry is known.
/// [Implicit](SuffixLen::Implicit) suffix lengths are recovered by extending
/// the k-mer as far as the encoder could have extended it at that position and
/// explicit and plain by taking the suffix length stored inside [SuffixLen].
///
/// Fails if the dictionary does not cover exactly `sequence_length` bases.
pub fn decode_path(
    path: Vec<(u32, SuffixLen)>,
    sequence_length: usize,
    sbwt: &SbwtIndexVariant,
) -> Result<Vec<u8>, GraphDecodeError> {
    let mut sequence: Vec<u8> = vec![0; sequence_length];
    let mut end = sequence_length;
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            for (colex_rank, suffix_len) in path {
                let kmer = sbwt.access_kmer(colex_rank as usize);
                let (kmer, len) = match suffix_len {
                    SuffixLen::Implicit => {
                        // The encoder extends at most end - k - 1 characters.
                        let new_kmer = crate::left_extend_kmer2(&kmer, sbwt, end.saturating_sub(k + 1));
                        let len = new_kmer.len();
                        (new_kmer, len)
                    },
                    SuffixLen::Explicit { suffix_len: len } => {
                        let new_kmer = crate::left_extend_kmer2(&kmer, sbwt, (len - k as u32) as usize);
                        assert_eq!(new_kmer.len(), len as usize);
                        (new_kmer, len as usize)
                    },
                    SuffixLen::Plain { suffix_len: len } => {
                        (kmer, len as usize)
                    },
                };
                if len > end {
                    return Err(GraphDecodeError(format!("dictionary covers more than {} bases", sequence_length)))
                }
                sequence[(end - len)..end].copy_from_slice(&kmer[(kmer.len() - len)..]);
                end -= len;
            }
        },
    };
    if end != 0 {
        return Err(GraphDecodeError(format!("dictionary covers {} of {} bases", sequence_length - end, sequence_length)))
    }

    Ok(sequence)
}

/// Builds the graph and remembers the path each colour takes.
pub struct GraphBuilder {

    /// colex rank -> node index (position in the set)
    nodes: IndexSet<u32>,

    /// (source node, target node, suffix length); edge id is the position in the set
    edges: IndexSet<(u32, u32, u32)>,

    /// Per edge id: true if the suffix length is implicit, see [SuffixLen]
    implicit: Vec<bool>,

    /// Per colour: edge ids (indexes into `edges`) in traversal order
    paths: Vec<Vec<usize>>,

    /// Store linear chains as single edges in [GraphBuilder::finish]
    collapse_chains: bool,
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
            edges: IndexSet::new(),
            implicit: Vec::new(),
            paths: Vec::new(),
            collapse_chains: false,
        }
    }

    /// Collapse linear chains into single edges in [GraphBuilder::finish].
    pub fn collapse_chains(
        mut self,
        collapse: bool,
    ) -> Self {
        self.collapse_chains = collapse;
        self
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
        implicit: bool,
    ) -> usize {
        let (id, inserted) = self.edges.insert_full((source, target, weight));
        if inserted {
            self.implicit.push(false);
        }
        // Implicit if any colour saw the extension stop in the SBWT: the length
        // then only depends on the _k_-mer, not on where the colour uses it.
        self.implicit[id] |= implicit;
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
        self.add_path_with_flags(sequence_length, dictionary, &vec![false; dictionary.len()])
    }

    /// Same as [GraphBuilder::add_path] but stores the suffix length of the
    /// entries flagged in `implicit` as [SuffixLen::Implicit].
    ///
    /// `dictionary` and `implicit` are the output of [crate::encode_sequence_flagged].
    pub fn add_path_with_flags(
        &mut self,
        sequence_length: usize,
        dictionary: &[(usize, Range<usize>)],
        implicit: &[bool],
    ) -> u32 {
        assert_eq!(dictionary.len(), implicit.len());
        let color: u32 = self.paths.len().try_into().unwrap();
        let mut path = Vec::with_capacity(dictionary.len() + 1);

        let mut prev = SOURCE_NODE;
        let mut prev_weight: u32 = sequence_length.try_into().expect("sequence length does not fit in u32");
        let mut prev_implicit = false;
        for ((suffix_len, colex), &is_implicit) in dictionary.iter().zip(implicit) {
            let node = self.insert_node(colex.start);
            path.push(self.edge(prev, node, prev_weight, prev_implicit));
            prev = node;
            prev_weight = (*suffix_len).try_into().unwrap();
            assert!(prev_weight != 0, "empty dictionary entry");
            prev_implicit = is_implicit;
        }
        path.push(self.edge(prev, SINK_NODE, prev_weight, prev_implicit));

        self.paths.push(path);
        color
    }

    /// Converts into the stored representation (graph + choice stream),
    /// collapsing linear chains if [GraphBuilder::collapse_chains] was set.
    pub fn finish(self) -> (StoredGraph, Choices) {
        if self.collapse_chains {
            let (chains, choices) = self.to_chains();
            (StoredGraph::Chains(chains), choices)
        } else {
            let (csr, choices) = self.to_csr();
            (StoredGraph::Csr(csr), choices)
        }
    }

    /// Converts into the stored representation (graph + choice stream).
    pub fn to_csr(self) -> (Csr<u32>, Choices) {
        let (csr, csr_pos) = self.build_csr();
        let choices = self.choices(&csr, &csr_pos);
        (csr, choices)
    }

    /// Converts into the stored representation with linear chains collapsed.
    ///
    /// Only junctions have more than one outgoing edge and their chains keep
    /// the CSR order of their first edges, so the choice stream is the same as
    /// for [GraphBuilder::to_csr].
    pub fn to_chains(self) -> (ChainCsr, Choices) {
        let (csr, csr_pos) = self.build_csr();
        let choices = self.choices(&csr, &csr_pos);
        (ChainCsr::from_csr(&csr), choices)
    }

    /// Records which edge each colour takes at every node with more than one outgoing edge.
    fn choices(
        &self,
        csr: &Csr<u32>,
        csr_pos: &[u32],
    ) -> Choices {
        let streams = self.paths.iter().map(|path| {
            let mut writer = ChoiceWriter::new();
            for &id in path {
                let source = self.edges.get_index(id).unwrap().0 as usize;
                let (start, end) = csr.row(source);
                if end - start > 1 {
                    writer.push((csr_pos[id] as usize - start) as u32, (end - start) as u32);
                }
            }
            writer.finish()
        }).collect();
        Choices { streams }
    }

    /// Builds the CSR and the CSR position of every edge id.
    fn build_csr(&self) -> (Csr<u32>, Vec<u32>) {
        let n_nodes = self.nodes.len();
        let n_edges = self.edges.len();

        // Number of times the colours take each edge.
        let mut uses = vec![0_u32; n_edges];
        for path in &self.paths {
            for &id in path {
                uses[id] += 1;
            }
        }

        // Edge order: by source, then the most used edges first, then target
        // and suffix length. The choice streams store ranks in this order, and
        // rank 0 being the most common choice makes them compress well.
        let mut order: Vec<usize> = (0..n_edges).collect();
        order.sort_unstable_by_key(|&id| {
            let &(source, target, weight) = self.edges.get_index(id).unwrap();
            (source, std::cmp::Reverse(uses[id]), target, weight)
        });

        let mut csr_pos = vec![0_u32; n_edges];
        let mut row_ptrs = vec![0_u32; n_nodes + 1];
        let mut column_indices = Vec::with_capacity(n_edges);
        let mut suffix_lens = Vec::with_capacity(n_edges);

        for (pos, &id) in order.iter().enumerate() {
            let &(source, target, weight) = self.edges.get_index(id).unwrap();
            csr_pos[id] = pos as u32;
            row_ptrs[source as usize + 1] += 1;
            column_indices.push(target);
            suffix_lens.push(if self.implicit[id] { IMPLICIT_SUFFIX_LEN } else { weight });
        }
        for i in 0..n_nodes {
            row_ptrs[i + 1] += row_ptrs[i];
        }

        let csr = Csr {
            node_weights: self.nodes.iter().copied().collect(),
            row_ptrs,
            column_indices,
            suffix_lens,
        };

        (csr, csr_pos)
    }
}

/// Picks the edge out of `node`, whose outgoing edges are `start..end`,
/// reading the choice stream if there is more than one.
fn take_edge(
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

/// Stored graph, either as is or with linear chains collapsed.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub enum StoredGraph {
    Csr(Csr<u32>),
    Chains(ChainCsr),
}

impl StoredGraph {
    /// Uncollapsed graph.
    pub fn into_csr(self) -> Csr<u32> {
        match self {
            StoredGraph::Csr(csr) => csr,
            StoredGraph::Chains(chains) => chains.to_csr(),
        }
    }

    /// One-line description for logging.
    pub fn summary(&self) -> String {
        match self {
            StoredGraph::Csr(csr) => format!(
                "{} nodes, {} edges",
                csr.node_weights.len(), csr.column_indices.len()),
            StoredGraph::Chains(chains) => format!(
                "{} junctions, {} chains, {} inner nodes, {} edges",
                chains.junction_weights.len(), chains.n_chains(), chains.inner_weights.len(), chains.suffix_lens.len()),
        }
    }
}

impl Walkable for StoredGraph {
    /// Walks the path of `color` from the source to the sink, see [walk_path].
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError> {
        match self {
            StoredGraph::Csr(csr) => csr.walk_path(choices, color),
            StoredGraph::Chains(chains) => chains.walk_path(choices, color),
        }
    }
}

/// Graph in compressed sparse row format; this is what gets stored.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct Csr<N> {
    pub node_weights: Vec<N>,
    pub row_ptrs: Vec<u32>,
    pub column_indices: Vec<u32>,
    pub suffix_lens: Vec<u32>,
}

impl<N> Csr<N> {
    /// Range of edge positions leaving `node`.
    pub fn row(
        &self,
        node: usize,
    ) -> (usize, usize) {
        (self.row_ptrs[node] as usize, self.row_ptrs[node + 1] as usize)
    }

    /// Walks the path of `color` from the source to the sink, calling
    /// `visit(node, edge position)` for every edge taken.
    ///
    /// The walk is bounded by the contig length, see [Walkable].
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

    /// Converts to a petgraph graph with `colorsets` from [Csr::colorsets] as the edge colours.
    pub fn to_petgraph<Ty: EdgeType>(
        self,
        colorsets: &[Vec<u32>],
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
            for (pos, colorset) in colorsets.iter().enumerate().take(end).skip(start) {
                graph.add_edge(
                    NodeIndex::new(source_idx),
                    NodeIndex::new(self.column_indices[pos] as usize),
                    ColexGraphEdge {
                        weight: self.suffix_lens[pos],
                        colors: colorset.clone(),
                    },
                );
            }
        }

        graph
    }
}

impl<N: Copy> Walkable for Csr<N> where u32: From<N> {
    /// Walks the path of `color` from the source to the sink.
    ///
    /// Returns the sequence length stored on the source edge and the NLZ
    /// dictionary as (colex rank, suffix length) pairs, in the order produced by
    /// [crate::encode_sequence].
    ///
    /// The returned value can be passed to [decode_path] to recover the nucleotide
    /// sequence.
    ///
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError> {
        let mut sequence_length = 0_u32;
        let mut dictionary: Vec<(u32, u32)> = Vec::new();

        self.walk(choices, color, |node, pos| {
            if node == SOURCE_NODE as usize {
                sequence_length = self.suffix_lens[pos];
            } else {
                let colex_rank: u32 = self.node_weights[node].into();
                dictionary.push((colex_rank, self.suffix_lens[pos]));
            }
        })?;

        Ok((sequence_length, dictionary))
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
