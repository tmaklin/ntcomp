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
//! ## Modules
//!
//! - [builder]: builds the graph from the contigs ([GraphBuilder]), with
//!   strand choice and edits;
//! - [csr]: the plain layout ([Csr]); [collapsed]: the chain layout and its
//!   stored sections;
//! - [paths]: walking a colour's path with its choices ([GraphWalk]);
//! - [choices]: the choice streams and the model that codes them;
//! - [snps]: the archive's SNP section;
//! - [encoder], [decoder]: the archive file; [header]: its header.
//!
pub mod builder;
pub mod choices;
pub mod collapsed;
pub mod csr;
pub mod decoder;
pub mod encoder;
pub mod header;
pub mod paths;
pub mod snps;

pub use builder::GraphBuilder;
pub use choices::Choices;
pub use csr::Csr;
pub use paths::GraphWalk;
use paths::take_edge;
use crate::walk::WalkIndex;
use collapsed::ChainCsr;

type E = Box<dyn std::error::Error>;

/// Bytes of a contig's blake3 hash kept in archives (enough to check a
/// decoded contig; the full 32 bytes were 6% of an archive).
pub const CONTIG_HASH_BYTES: usize = 8;

/// A contig's hash as stored in archives: the first [CONTIG_HASH_BYTES] bytes
/// of its blake3 hash.
pub type ContigHash = [u8; CONTIG_HASH_BYTES];

/// [ContigHash] of a (normalized) contig sequence.
pub fn contig_hash(
    sequence: &[u8],
) -> ContigHash {
    let mut out = [0u8; CONTIG_HASH_BYTES];
    out.copy_from_slice(&blake3::hash(sequence).as_bytes()[..CONTIG_HASH_BYTES]);
    out
}

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

/// Serialized form of [StoredGraph] (the same bytes as its serde form).
#[derive(serde::Serialize)]
enum StoredGraphRef<'a> {
    Csr(&'a Csr<u32>),
    Chains(collapsed::StoredChains),
}

/// Deserialized form of [StoredGraph] before the chains are decoded.
#[derive(serde::Deserialize)]
enum StoredGraphOwned {
    Csr(Csr<u32>),
    Chains(collapsed::StoredChains),
}

/// Stored graph, either as is or with linear chains collapsed.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub enum StoredGraph {
    Csr(Csr<u32>),
    Chains(ChainCsr),
}

impl StoredGraph {
    /// Serializes the graph. With the `index` (and its LCS array), the
    /// suffix lengths of collapsed chains are coded per LCS context and their
    /// nodes are derived by walks in the index where possible, see
    /// [collapsed]; reading them back then needs the index. `k` is the
    /// _k_-mer size of the index (None if unknown), which sets the direct
    /// suffix length symbols, see [collapsed::max_direct_suffix_len].
    pub fn to_bytes(
        &self,
        index: Option<WalkIndex>,
        k: Option<u32>,
    ) -> Result<Vec<u8>, postcard::Error> {
        let repr = match self {
            StoredGraph::Csr(csr) => StoredGraphRef::Csr(csr),
            StoredGraph::Chains(chains) => StoredGraphRef::Chains(chains.to_stored_with(index, collapsed::max_direct_suffix_len(k), k, index.is_some())),
        };
        postcard::to_allocvec(&repr)
    }

    /// Reads a graph written by [StoredGraph::to_bytes]; `index` is needed if
    /// it was written with the index (its LCS array is optional: the values
    /// are otherwise computed from the SBWT, slower), and `k` must be the
    /// _k_-mer size it was written with (None works for _k_ ≤ 255).
    pub fn from_bytes(
        bytes: &[u8],
        index: Option<WalkIndex>,
        k: Option<u32>,
    ) -> Result<Self, E> {
        Ok(match postcard::from_bytes(bytes)? {
            StoredGraphOwned::Csr(csr) => StoredGraph::Csr(csr),
            StoredGraphOwned::Chains(stored) => StoredGraph::Chains(ChainCsr::from_stored_with(stored, index, k)?),
        })
    }

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

impl GraphWalk for StoredGraph {
    fn n_nodes(&self) -> usize {
        match self { StoredGraph::Csr(g) => g.n_nodes(), StoredGraph::Chains(g) => g.n_nodes() }
    }

    fn n_edges(&self) -> usize {
        match self { StoredGraph::Csr(g) => g.n_edges(), StoredGraph::Chains(g) => g.n_edges() }
    }

    fn row(&self, node: usize) -> (usize, usize) {
        match self { StoredGraph::Csr(g) => g.row(node), StoredGraph::Chains(g) => g.row(node) }
    }

    fn target(&self, edge: usize) -> usize {
        match self { StoredGraph::Csr(g) => g.target(edge), StoredGraph::Chains(g) => g.target(edge) }
    }

    fn max_steps(&self, edge: usize) -> usize {
        match self { StoredGraph::Csr(g) => g.max_steps(edge), StoredGraph::Chains(g) => g.max_steps(edge) }
    }

    /// Walks the path of `color` from the source to the sink, see [GraphWalk::walk_path].
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

#[derive(Debug)]
pub struct GraphDecodeError(pub String);

impl std::fmt::Display for GraphDecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "graph decoding failed: {}", self.0)
    }
}

impl std::error::Error for GraphDecodeError {}
