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
//! - [encoder], [decoder]: the archive file.
//!
use sbwt::LcsArray;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

pub mod builder;
pub mod choices;
pub mod collapsed;
pub mod csr;
pub mod decoder;
pub mod encoder;
pub mod paths;
pub mod snps;

pub use builder::GraphBuilder;
pub use choices::Choices;
pub use csr::Csr;
pub use paths::GraphWalk;
use paths::take_edge;
use crate::walk::IndexWalk;
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

/// [Implicit](SuffixLen::Implicit) suffix lengths are stored as a 0; a real suffix length can never be 0.
const IMPLICIT_SUFFIX_LEN: u32 = 0;

/// Suffix length and type for an edge
///
/// Classifies the suffix length for a dictionary entry by specifying whether
/// its implicit (for a left extended _k_-mer but can be inferred), explicit
/// (left extended and needs to be declared), or plain (subsequence of a plan
/// _k_-mer). Stores the value for explicit and plain types.
///
pub enum SuffixLen {
    /// The dictionary entry is a left extended _k_-mer whose extension stopped
    /// where [crate::left_extend_kmer] stops, so [decode_path] recovers the
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
/// dictionary (colex rank, suffix length), as returned by [GraphWalk::walk_path].
///
/// The dictionary is processed in its own order, from the end of the sequence
/// to the start, so that the end position of every entry is known.
/// [Implicit](SuffixLen::Implicit) suffix lengths are recovered by extending
/// the k-mer as far as the encoder could have extended it at that position and
/// explicit and plain by taking the suffix length stored inside [SuffixLen].
///
/// The k-mers are extended by walking the de Bruijn graph of the SBWT
/// ([IndexWalk::left_extend](crate::walk::IndexWalk::left_extend)), with
/// `lcs` if given and with LCS values computed from the SBWT otherwise.
///
/// Fails if the dictionary does not cover exactly `sequence_length` bases.
pub fn decode_path(
    path: Vec<(u32, SuffixLen)>,
    sequence_length: usize,
    sbwt: &SbwtIndexVariant,
    lcs: Option<&LcsArray>,
) -> Result<Vec<u8>, GraphDecodeError> {
    let mut sequence: Vec<u8> = vec![0; sequence_length];
    let mut end = sequence_length;
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            // Left extension of the k-mer at `colex_rank` by walking the SBWT.
            let walk = crate::walk::WalkIndex { sbwt, lcs };
            let extend = |colex_rank: u32, max_extension_len: usize, implicit: bool| walk.left_extend(colex_rank as usize, max_extension_len, implicit).0;
            for (colex_rank, suffix_len) in path {
                let (kmer, len) = match suffix_len {
                    SuffixLen::Implicit => {
                        // The encoder extends at most end - k - 1 characters.
                        let new_kmer = extend(colex_rank, end.saturating_sub(k + 1), true);
                        let len = new_kmer.len();
                        (new_kmer, len)
                    },
                    SuffixLen::Explicit { suffix_len: len } => {
                        let new_kmer = extend(colex_rank, (len - k as u32) as usize, false);
                        assert_eq!(new_kmer.len(), len as usize);
                        (new_kmer, len as usize)
                    },
                    SuffixLen::Plain { suffix_len: len } => {
                        (sbwt.access_kmer(colex_rank as usize), len as usize)
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

/// Source of LCS array values for decoding chains stored with LCS contexts,
/// see [collapsed].
#[derive(Clone, Copy)]
pub enum LcsSource<'a> {
    /// The LCS array of the index.
    Array(&'a LcsArray),

    /// Computed from the index with [crate::walk::lcs_from_sbwt] (slower).
    Sbwt(&'a SbwtIndexVariant),
}

impl LcsSource<'_> {
    /// LCS array value of `colex_rank`.
    pub fn lcs(
        &self,
        colex_rank: u32,
    ) -> u32 {
        match self {
            LcsSource::Array(lcs) => lcs.access(colex_rank as usize) as u32,
            LcsSource::Sbwt(SbwtIndexVariant::SubsetMatrix(sbwt)) => crate::walk::lcs_from_sbwt(sbwt, colex_rank as usize) as u32,
        }
    }
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
    /// Serializes the graph. With `lcs`, the suffix lengths of collapsed
    /// chains are coded per LCS context, and with `sbwt` too their nodes are
    /// derived by walks in the index where possible, see [collapsed]; reading
    /// them back then needs an [LcsSource] (and the SBWT). `k` is the _k_-mer
    /// size of the index (None if unknown), which sets the direct suffix
    /// length symbols, see [collapsed::max_direct_suffix_len].
    pub fn to_bytes(
        &self,
        lcs: Option<&LcsArray>,
        sbwt: Option<&SbwtIndexVariant>,
        k: Option<u32>,
    ) -> Result<Vec<u8>, postcard::Error> {
        let repr = match self {
            StoredGraph::Csr(csr) => StoredGraphRef::Csr(csr),
            StoredGraph::Chains(chains) => match lcs {
                Some(lcs) => {
                    let lookup = |rank: u32| lcs.access(rank as usize) as u32;
                    let walk = sbwt.map(|SbwtIndexVariant::SubsetMatrix(sbwt)| collapsed::WalkIndex { sbwt, lcs: Some(lcs) });
                    StoredGraphRef::Chains(chains.to_stored_with(Some(&lookup), collapsed::max_direct_suffix_len(k), k, walk))
                },
                None => StoredGraphRef::Chains(chains.to_stored_with(None, 0, k, None)),
            },
        };
        postcard::to_allocvec(&repr)
    }

    /// Reads a graph written by [StoredGraph::to_bytes]; `lcs` is needed if
    /// it was written with the LCS array, `sbwt` if it was written with the
    /// SBWT too, and `k` must be the _k_-mer size it was written with (None
    /// works for _k_ ≤ 255).
    pub fn from_bytes(
        bytes: &[u8],
        lcs: Option<LcsSource>,
        sbwt: Option<&SbwtIndexVariant>,
        k: Option<u32>,
    ) -> Result<Self, E> {
        Ok(match postcard::from_bytes(bytes)? {
            StoredGraphOwned::Csr(csr) => StoredGraph::Csr(csr),
            StoredGraphOwned::Chains(stored) => {
                let lookup = lcs.map(|source| move |rank: u32| source.lcs(rank));
                let lookup: Option<collapsed::LcsLookup> = lookup.as_ref().map(|f| f as collapsed::LcsLookup);
                let array = match lcs { Some(LcsSource::Array(array)) => Some(array), _ => None };
                let walk = sbwt.map(|SbwtIndexVariant::SubsetMatrix(sbwt)| collapsed::WalkIndex { sbwt, lcs: array });
                StoredGraph::Chains(ChainCsr::from_stored_with(stored, lookup, k, walk)?)
            },
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
