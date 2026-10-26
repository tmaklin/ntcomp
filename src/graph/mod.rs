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
//! Linear chains of nodes are stored as single edges between junctions
//! ([ChainCsr]), see [collapsed]. Linear chains are paths on the graph composed
//! entirely of inner nodes, meaning nodes that have exactly one incoming and
//! exactly one outgoing edge. All other nodes are junction nodes.
//!
//! ## Modules
//!
//! - [builder]: builds the graph from the contigs ([GraphBuilder]), with
//!   strand choice and edits;
//! - [csr]: the uncollapsed graph ([Csr]), built first and expanded for
//!   viewing; [collapsed]: the stored chain layout and its sections;
//! - [paths]: walking a colour's path with its choices ([GraphWalk]);
//! - [choices]: the choice streams and the model that codes them;
//! - [snps]: the archive's SNP section; [names]: its contig names;
//! - [encoder], [decoder]: the archive file; [header]: its header.
//!
pub mod builder;
pub mod choices;
pub mod collapsed;
pub mod csr;
pub mod decoder;
pub mod encoder;
pub mod header;
pub mod names;
pub mod paths;
pub mod snps;

pub use builder::GraphBuilder;
pub use choices::Choices;
pub use collapsed::ChainCsr;
pub use csr::Csr;
pub use paths::GraphWalk;
use paths::take_edge;
use crate::nlz::Entry;
use crate::walk::WalkIndex;

type E = Box<dyn std::error::Error>;

/// Bytes of an input file's hash kept in archives (enough to check a
/// decoded file; per-contig hashes were a fifth of an archive of fragmented
/// assemblies).
pub const FILE_HASH_BYTES: usize = 8;

/// An input file's hash as stored in archives, see [file_hash].
pub type FileHash = [u8; FILE_HASH_BYTES];

/// A contig's blake3 hash, see [contig_digest].
pub type ContigDigest = [u8; 32];

/// The blake3 hash of a contig: the length of its name (u64,
/// little-endian), its name and its (normalized) sequence.
pub fn contig_digest(
    name: &[u8],
    sequence: &[u8],
) -> ContigDigest {
    let mut hasher = blake3::Hasher::new();
    hasher.update(&(name.len() as u64).to_le_bytes());
    hasher.update(name);
    hasher.update(sequence);
    *hasher.finalize().as_bytes()
}

/// The [FileHash] of a file's contigs: the first [FILE_HASH_BYTES] bytes of
/// the blake3 hash of their [contig_digest]s in order, so that the contigs
/// can be checked in any order.
pub fn file_hash(
    digests: &[ContigDigest],
) -> FileHash {
    let mut hasher = blake3::Hasher::new();
    for digest in digests {
        hasher.update(digest);
    }
    hasher.finalize().as_bytes()[..FILE_HASH_BYTES].try_into().unwrap()
}

/// Node index of the dummy source node.
const SOURCE_NODE: u32 = 0;
/// Node index of the dummy sink node.
const SINK_NODE: u32 = 1;
/// Colex rank stored in the dummy source node. Real colex ranks must be smaller.
const SOURCE_COLEX: u32 = u32::MAX - 1;
/// Colex rank stored in the dummy sink node. Real colex ranks must be smaller.
const SINK_COLEX: u32 = u32::MAX;

/// The length of an [implicit](Entry::implicit) dictionary entry as stored on
/// an edge (a real length is never 0).
const IMPLICIT_SUFFIX_LEN: u32 = 0;

/// The dictionary entry of the node `colex_rank` whose edge stores `len`.
fn stored_entry(
    colex_rank: u32,
    len: u32,
) -> Entry {
    Entry { colex_rank, len, implicit: len == IMPLICIT_SUFFIX_LEN }
}

/// Leading byte of the stored graph: the chain layout (the variant tag that
/// separated it from the plain layout, which is no longer written or read).
const CHAINS_TAG: u8 = 1;

impl ChainCsr {
    /// Serializes the graph. With the `index` (and its LCS array), the
    /// suffix lengths are coded per LCS context and the nodes are derived by
    /// walks in the index where possible, see [collapsed]; reading them back
    /// then needs the index. `k` is the _k_-mer size of the index (None if
    /// unknown), which sets the direct suffix length symbols, see
    /// [collapsed::max_direct_suffix_len].
    pub fn to_bytes(
        &self,
        index: Option<WalkIndex>,
        k: Option<u32>,
    ) -> Result<Vec<u8>, postcard::Error> {
        let stored = self.to_stored_with(index, collapsed::max_direct_suffix_len(k), k, index.is_some());
        let mut bytes = vec![CHAINS_TAG];
        bytes.extend(postcard::to_allocvec(&stored)?);
        Ok(bytes)
    }

    /// Reads a graph written by [ChainCsr::to_bytes]; `index` is needed if
    /// it was written with the index (its LCS array is optional: the values
    /// are otherwise computed from the SBWT, slower), and `k` must be the
    /// _k_-mer size it was written with (None works for _k_ ≤ 255).
    pub fn from_bytes(
        bytes: &[u8],
        index: Option<WalkIndex>,
        k: Option<u32>,
    ) -> Result<Self, E> {
        match bytes.split_first() {
            Some((&CHAINS_TAG, rest)) => Ok(ChainCsr::from_stored_with(postcard::from_bytes(rest)?, index, k)?),
            Some((0, _)) => Err(Box::new(GraphDecodeError("the plain graph layout is no longer supported".to_string()))),
            _ => Err(Box::new(GraphDecodeError("unknown graph layout".to_string()))),
        }
    }

    /// One-line description for logging.
    pub fn summary(&self) -> String {
        format!(
            "{} junctions, {} chains, {} inner nodes, {} edges",
            self.junction_weights.len(), self.n_chains(), self.inner_weights.len(), self.suffix_lens.len())
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
