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
//! NLZ-OG graph with linear chains collapsed into single edges.
//!
//! ## Graph structure
//!
//! Represents the graph as paths through inner nodes and junctions. An *inner
//! node* is a node with exactly one incoming and one outgoing edge. Every other
//! node including the dummy source and sink are *junction nodes*. Following
//! edges from a junction through inner nodes always ends at a junction, so the
//! graph can be stored as linear chains between junctions.
//!
//! Node weights for junction nodes are colex ranks. Node weights for the
//! collapsed linear chains between junctions are concatenated in the edge
//! weight.
//!
//! ## Colours
//!
//! All edges between inner nodes always have the same colour set, as changing the
//! colour set requires introducing an incoming edge with a different colour
//! set, which converts the node into a junction. Therefore, a chain stores the
//! colex ranks of its inner nodes and the suffix lengths of its edges, and its
//! colour set only once.
//!
//! Colours only have a choice at junctions, so the choice stream has the same
//! number of entries as for the uncollapsed graph. The stored ranks refer to
//! the order of the chains in the junction's row. This means that the edges
//! between inner nodes do not need to store the colour sets if they are stored
//! in the choice stream.
//!
//! ## Stored colex ranks
//!
//! The colex ranks of the junctions and the inner nodes are bit-packed with the
//! same width for every rank. The number of bits needed for the largest rank
//! depends on the size of the SBWT. The dummy source and sink are always
//! junctions 0 and 1 and are not stored.
//!
//! ## Stored suffix lengths
//!
//! The suffix lengths on the edges out of the dummy source node are the total
//! contig lengths and are stored as they are.
//!
//! Every other suffix length is a symbol of a canonical Huffman code, in
//! chain order:
//! - 0: the length is implicit, see [SuffixLen::Implicit](super::SuffixLen),
//! - 1..=D: the length itself,
//! - D + 1: a longer length, followed by the length - D (Elias gamma code),
//!
//! where D = max(_k_, [MIN_DIRECT_SUFFIX_LEN]) (see [max_direct_suffix_len]),
//! with _k_ taken from the SBWT when the graph is stored or read with the
//! index, so every short entry (< _k_) has its own symbol. Without the index D
//! is [MIN_DIRECT_SUFFIX_LEN], the same as for any _k_ ≤ 255.
//!
//! Without the LCS array there is one code for all edges. With it, there is
//! one code per context min(LCS(r), D), where r is the colex
//! rank of the edge's source node: a short entry's rank is the start of the
//! SBWT interval of its suffix, so its length is always larger than LCS(r),
//! often by one. The code lengths of all contexts are LZMA compressed. The
//! decoder takes LCS(r) from the LCS array or computes it from the SBWT with
//! [crate::lcs_from_sbwt].
//!
use super::{
    Csr,
    Choices,
    GraphDecodeError,
    SINK_COLEX,
    SINK_NODE,
    SOURCE_COLEX,
    SOURCE_NODE,
    take_edge,
    Walkable,
};
use super::bits::{BitReader, BitWriter};
use super::huffman::HuffmanCode;

/// Smallest number of suffix lengths with their own Huffman symbol: lengths
/// up to this are coded directly even when _k_ is smaller (left extensions
/// longer than _k_ code better directly than escaped).
pub const MIN_DIRECT_SUFFIX_LEN: u32 = 255;

/// Longest suffix length with its own Huffman symbol for _k_-mer size `k`
/// (None: unknown), see the module docs. Also the cap of the LCS contexts.
pub fn max_direct_suffix_len(
    k: Option<u32>,
) -> u32 {
    k.unwrap_or(0).max(MIN_DIRECT_SUFFIX_LEN)
}

/// LCS array value of a colex rank, see [super::LcsSource].
pub type LcsLookup<'a> = &'a (dyn Fn(u32) -> u32 + Sync);

/// Graph of linear chains between junctions, in compressed sparse row format.
///
/// Serialized as [StoredChains].
#[derive(Debug, PartialEq, Eq)]
pub struct ChainCsr {
    /// Colex ranks of the junction nodes.
    ///
    /// Junction with weight 0 is the dummy source node, weight 1 is the sink.
    pub junction_weights: Vec<u32>,

    /// Chains leaving junction node `j` are `row_ptrs[j]..row_ptrs[j + 1]`.
    pub row_ptrs: Vec<u32>,

    /// Index of the junction that each chain ends at.
    ///
    /// Junctions are numbered in breadth-first order.
    /// TODO describe breadth-first order in the module docs.
    pub targets: Vec<u32>,

    /// Inner nodes of linear chain `c` are
    /// `inner_weights[chain_ptrs[c]..chain_ptrs[c + 1]]`.
    pub chain_ptrs: Vec<u32>,

    /// Colex ranks of the inner nodes, concatenated over all chains.
    ///
    /// Colex ranks for chain `c` are stored at
    /// `inner_weights[chain_ptrs[c]..chain_ptrs[c + 1]]`.
    pub inner_weights: Vec<u32>,

    /// Suffix lengths of the edges in each chain, concatenated over all chains.
    ///
    /// A chain with `n` inner nodes has `n + 1` edges, so chain `c` starts at
    /// `chain_ptrs[c] + c`.
    pub suffix_lens: Vec<u32>,
}

/// Stored form of [ChainCsr], see the module docs.
#[derive(serde::Serialize, serde::Deserialize, Debug)]
pub struct StoredChains {
    /// Number of junction nodes, including the dummy source and sink.
    pub n_junctions: u32,

    /// Number of bits used for every stored colex rank.
    pub colex_bits: u8,

    /// Colex ranks of the junction nodes except the dummy source and sink,
    /// `colex_bits` bits each.
    pub junction_weights: Vec<u8>,

    /// Colex ranks of the inner nodes, concatenated over all chains,
    /// `colex_bits` bits each.
    pub inner_weights: Vec<u8>,

    /// Suffix lengths of the edges out of the dummy source node (the contig lengths).
    pub sequence_lens: Vec<u32>,

    /// Largest LCS context of the suffix length codes; 0: one code for all.
    ///
    /// Otherwise there is one code per context min(LCS(source rank), lcs_cap)
    /// and `suffix_code` holds their code lengths, LZMA compressed.
    pub lcs_cap: u16,

    /// Huffman code lengths of the other suffix length symbols.
    pub suffix_code: Vec<u8>,

    /// Other suffix lengths, Huffman coded.
    pub suffix_lens: Vec<u8>,

    /// Out-degree of every junction except the sink.
    pub out_degrees: Vec<u8>,

    /// Number of inner nodes + 1 of every chain.
    pub chain_lens: Vec<u8>,

    /// Target of every chain.
    pub targets: Vec<u8>,
}

/// Number of bits needed for colex ranks up to `max_colex` (at least 1).
fn colex_bits(
    max_colex: u32,
) -> u32 {
    (32 - max_colex.leading_zeros()).max(1)
}

/// Packs `ranks` with `bits` bits each.
fn pack_colex(
    ranks: &[u32],
    bits: u32,
) -> Vec<u8> {
    let mut writer = BitWriter::new();
    for &rank in ranks {
        writer.write_bits(rank as u64, bits);
    }
    writer.finish()
}

/// Unpacks `n` ranks of `bits` bits each; None unless `bytes` holds exactly that many.
fn unpack_colex(
    bytes: &[u8],
    n: usize,
    bits: u32,
) -> Option<Vec<u32>> {
    if bytes.len() != (n * bits as usize).div_ceil(8) {
        return None
    }
    let mut reader = BitReader::new(bytes);
    (0..n).map(|_| reader.read_bits(bits).map(|rank| rank as u32)).collect()
}

impl ChainCsr {
    /// Collapses the linear chains in `csr`.
    pub fn from_csr(
        csr: &Csr<u32>,
    ) -> Self {
        let n_nodes = csr.node_weights.len();
        let n_edges = csr.column_indices.len();

        let mut in_degree = vec![0_u32; n_nodes];
        for &target in &csr.column_indices {
            in_degree[target as usize] += 1;
        }

        let is_junction: Vec<bool> = (0..n_nodes).map(|node| {
            let (start, end) = csr.row(node);
            let is_dummy = node == SOURCE_NODE as usize || node == SINK_NODE as usize;
            is_dummy || in_degree[node] != 1 || end - start != 1
        }).collect();

        // Node index -> junction index in breadth-first order, u32::MAX until
        // reached. `order` is the queue of the search and lists the nodes of
        // the junctions in index order.
        let mut junction_index = vec![u32::MAX; n_nodes];
        junction_index[SOURCE_NODE as usize] = 0;
        junction_index[SINK_NODE as usize] = 1;
        let mut order: Vec<usize> = vec![SOURCE_NODE as usize, SINK_NODE as usize];

        let mut row_ptrs = Vec::with_capacity(n_nodes + 1);
        let mut targets = Vec::new();
        let mut chain_ptrs = vec![0_u32];
        let mut inner_weights = Vec::new();
        let mut suffix_lens = Vec::with_capacity(n_edges);

        row_ptrs.push(0);
        let mut head = 0;
        while head < order.len() {
            let node = order[head];
            head += 1;

            let (start, end) = csr.row(node);
            for pos in start..end {
                suffix_lens.push(csr.suffix_lens[pos]);

                let mut next = csr.column_indices[pos] as usize;
                while !is_junction[next] {
                    let edge = csr.row(next).0;
                    inner_weights.push(csr.node_weights[next]);
                    suffix_lens.push(csr.suffix_lens[edge]);
                    next = csr.column_indices[edge] as usize;
                }

                if junction_index[next] == u32::MAX {
                    junction_index[next] = order.len() as u32;
                    order.push(next);
                }
                targets.push(junction_index[next]);
                chain_ptrs.push(inner_weights.len() as u32);
            }
            row_ptrs.push(targets.len() as u32);
        }
        assert_eq!(order.len(), is_junction.iter().filter(|x| **x).count(), "junction not reachable from the source");

        ChainCsr {
            junction_weights: order.iter().map(|&node| csr.node_weights[node]).collect(),
            row_ptrs,
            targets,
            chain_ptrs,
            inner_weights,
            suffix_lens,
        }
    }

    /// Converts into the stored form.
    ///
    /// Panics if the junctions are not numbered in breadth-first order.
    pub fn to_stored(
        &self,
    ) -> StoredChains {
        self.to_stored_with(None, 0, None)
    }

    /// Colex rank of the source of every edge, in chain order (the junction
    /// for the first edge of a chain, else the previous inner node).
    fn edge_sources(
        &self,
    ) -> Vec<u32> {
        let mut sources = Vec::with_capacity(self.suffix_lens.len());
        for junction in 0..self.junction_weights.len() {
            let (start, end) = self.row(junction);
            for chain in start..end {
                sources.push(self.junction_weights[junction]);
                sources.extend_from_slice(self.inner_nodes(chain));
            }
        }
        sources
    }

    /// Same as [ChainCsr::to_stored], with one suffix length code per LCS
    /// value of the edge's source rank (capped at `lcs_cap`) if `lcs` gives
    /// the LCS array value of a colex rank, and direct suffix length symbols
    /// for _k_-mer size `k`, see the module docs.
    pub fn to_stored_with(
        &self,
        lcs: Option<LcsLookup>,
        lcs_cap: u32,
        k: Option<u32>,
    ) -> StoredChains {
        let lcs_cap = if lcs.is_some() { lcs_cap.max(1) } else { 0 };
        let max_direct = max_direct_suffix_len(k);
        let long_symbol = max_direct as usize + 1;
        let n_junctions = self.junction_weights.len();

        let mut out_degrees = BitWriter::new();
        for junction in (0..n_junctions).filter(|&j| j != SINK_NODE as usize) {
            let (start, end) = self.row(junction);
            out_degrees.write_gamma((end - start) as u64);
        }
        let (sink_start, sink_end) = self.row(SINK_NODE as usize);
        assert_eq!(sink_start, sink_end, "sink has outgoing chains");

        let mut chain_lens = BitWriter::new();
        let mut targets = BitWriter::new();
        let mut next_new = 2_u32;
        for chain in 0..self.n_chains() {
            chain_lens.write_gamma(self.inner_nodes(chain).len() as u64 + 1);

            let target = self.targets[chain];
            if target == next_new {
                targets.write_bits(1, 1);
                next_new += 1;
            } else {
                assert!(target < next_new, "junctions are not numbered in breadth-first order");
                targets.write_bits(0, 1);
                targets.write_truncated_binary(target, next_new);
            }
        }
        assert_eq!(next_new as usize, n_junctions, "junctions are not numbered in breadth-first order");

        // Chains out of the source come first; their first edge holds the contig length.
        let n_source_chains = self.row(SOURCE_NODE as usize).1;
        let sources = if lcs_cap > 0 { self.edge_sources() } else { Vec::new() };
        let mut sequence_lens = Vec::with_capacity(n_source_chains);
        // (context, symbol, length)
        let mut symbols: Vec<(usize, usize, u32)> = Vec::with_capacity(self.suffix_lens.len());
        let mut edge = 0;
        for chain in 0..self.n_chains() {
            for (i, &len) in self.chain_suffix_lens(chain).iter().enumerate() {
                let source = if lcs_cap > 0 { sources[edge] } else { 0 };
                edge += 1;
                if chain < n_source_chains && i == 0 {
                    sequence_lens.push(len);
                    continue;
                }
                let context = if lcs_cap > 0 { lcs.unwrap()(source).min(lcs_cap) as usize } else { 0 };
                if len > max_direct {
                    symbols.push((context, long_symbol, len));
                } else {
                    symbols.push((context, len as usize, len));
                }
            }
        }
        let n_contexts = lcs_cap as usize + 1;
        let mut counts = vec![vec![0_u64; long_symbol + 1]; n_contexts];
        for &(context, symbol, _) in &symbols {
            counts[context][symbol] += 1;
        }
        let codes: Vec<HuffmanCode> = counts.iter().map(|c| HuffmanCode::from_counts(c)).collect();
        let mut suffix_lens = BitWriter::new();
        for &(context, symbol, len) in &symbols {
            codes[context].write(symbol, &mut suffix_lens);
            if symbol == long_symbol {
                suffix_lens.write_gamma((len - max_direct) as u64);
            }
        }
        let suffix_code = if lcs_cap == 0 {
            codes[0].lengths().to_vec()
        } else {
            let lengths: Vec<u8> = codes.iter().flat_map(|c| c.lengths().iter().copied()).collect();
            postcard::to_allocvec(&super::lzma_wrapper::compress(&lengths).unwrap()).unwrap()
        };

        // Real colex ranks only; the dummy source and sink are junctions 0 and 1.
        assert_eq!(self.junction_weights[..2], [SOURCE_COLEX, SINK_COLEX], "junctions 0 and 1 are not the source and sink");
        let real_junctions = &self.junction_weights[2..];
        let max_colex = real_junctions.iter().chain(&self.inner_weights).copied().max().unwrap_or(0);
        assert!(max_colex < SOURCE_COLEX, "colex rank collides with dummy node ranks");
        let colex_bits = colex_bits(max_colex);

        StoredChains {
            n_junctions: n_junctions.try_into().expect("too many junctions"),
            colex_bits: colex_bits as u8,
            junction_weights: pack_colex(real_junctions, colex_bits),
            inner_weights: pack_colex(&self.inner_weights, colex_bits),
            sequence_lens,
            lcs_cap: lcs_cap as u16,
            suffix_code,
            suffix_lens: suffix_lens.finish(),
            out_degrees: out_degrees.finish(),
            chain_lens: chain_lens.finish(),
            targets: targets.finish(),
        }
    }

    /// Converts from the stored form.
    pub fn from_stored(
        stored: StoredChains,
    ) -> Result<Self, GraphDecodeError> {
        Self::from_stored_with(stored, None, None)
    }

    /// Same as [ChainCsr::from_stored]; `lcs` must give the LCS array value
    /// of a colex rank if the chains were stored with LCS contexts, and `k`
    /// must be the _k_-mer size they were stored with (None for _k_ ≤ 255).
    pub fn from_stored_with(
        stored: StoredChains,
        lcs: Option<LcsLookup>,
        k: Option<u32>,
    ) -> Result<Self, GraphDecodeError> {
        let max_direct = max_direct_suffix_len(k);
        let long_symbol = max_direct as usize + 1;
        let truncated = |what: &str| GraphDecodeError(format!("stored chains: {} truncated or invalid", what));
        let n_junctions = stored.n_junctions as usize;
        if n_junctions < 2 {
            return Err(GraphDecodeError("stored chains: missing source or sink".to_string()))
        }
        let colex_bits = stored.colex_bits as u32;
        if !(1..=32).contains(&colex_bits) {
            return Err(truncated("colex rank width"))
        }
        let mut junction_weights = vec![SOURCE_COLEX, SINK_COLEX];
        junction_weights.extend(unpack_colex(&stored.junction_weights, n_junctions - 2, colex_bits).ok_or_else(|| truncated("junction colex ranks"))?);

        let mut row_ptrs = Vec::with_capacity(n_junctions + 1);
        row_ptrs.push(0_u32);
        let mut degrees = BitReader::new(&stored.out_degrees);
        for junction in 0..n_junctions {
            let degree = if junction == SINK_NODE as usize {
                0
            } else {
                degrees.read_gamma().ok_or_else(|| truncated("out-degrees"))?
            };
            let ptr = (*row_ptrs.last().unwrap() as u64 + degree).try_into().map_err(|_| truncated("out-degrees"))?;
            row_ptrs.push(ptr);
        }
        let n_chains = *row_ptrs.last().unwrap() as usize;

        let mut chain_ptrs = Vec::with_capacity(n_chains + 1);
        chain_ptrs.push(0_u32);
        let mut lens = BitReader::new(&stored.chain_lens);
        for _ in 0..n_chains {
            let len = lens.read_gamma().ok_or_else(|| truncated("chain lengths"))? - 1;
            let ptr = (*chain_ptrs.last().unwrap() as u64 + len).try_into().map_err(|_| truncated("chain lengths"))?;
            chain_ptrs.push(ptr);
        }
        let n_inner = *chain_ptrs.last().unwrap() as usize;
        let inner_weights = unpack_colex(&stored.inner_weights, n_inner, colex_bits).ok_or_else(|| truncated("inner colex ranks"))?;

        let n_source_chains = row_ptrs[1] as usize;
        if stored.sequence_lens.len() != n_source_chains {
            return Err(truncated("contig lengths"))
        }
        let lcs_cap = stored.lcs_cap as u32;
        let codes: Vec<HuffmanCode> = if lcs_cap == 0 {
            if stored.suffix_code.len() != long_symbol + 1 {
                return Err(GraphDecodeError(format!("stored chains: suffix code for {} symbols, expected {} (stored with another k?)", stored.suffix_code.len(), long_symbol + 1)))
            }
            vec![HuffmanCode::from_lengths(stored.suffix_code).map_err(|e| GraphDecodeError(format!("stored chains: {}", e)))?]
        } else {
            let lzma: super::lzma_wrapper::LzmaData = postcard::from_bytes(&stored.suffix_code).map_err(|_| truncated("suffix codes"))?;
            let lengths = super::lzma_wrapper::decompress(&lzma).map_err(|_| truncated("suffix codes"))?;
            if lengths.len() != (lcs_cap as usize + 1) * (long_symbol + 1) {
                return Err(truncated("suffix codes"))
            }
            lengths.chunks(long_symbol + 1).map(|l| HuffmanCode::from_lengths(l.to_vec())).collect::<Result<_, _>>().map_err(|e| GraphDecodeError(format!("stored chains: {}", e)))?
        };
        let lcs = match (lcs_cap, lcs) {
            (0, _) => None,
            (_, Some(lcs)) => Some(lcs),
            (_, None) => return Err(GraphDecodeError("stored chains: the suffix lengths need the LCS array".to_string())),
        };
        let mut bits = BitReader::new(&stored.suffix_lens);
        let mut suffix_lens = Vec::with_capacity(n_inner + n_chains);
        let mut junction = 0;
        for chain in 0..n_chains {
            while row_ptrs[junction + 1] as usize <= chain {
                junction += 1;
            }
            let n_edges = (chain_ptrs[chain + 1] - chain_ptrs[chain]) as usize + 1;
            for i in 0..n_edges {
                if chain < n_source_chains && i == 0 {
                    suffix_lens.push(stored.sequence_lens[chain]);
                    continue;
                }
                let code = match lcs {
                    None => &codes[0],
                    Some(lcs) => {
                        let source = if i == 0 { junction_weights[junction] } else { inner_weights[chain_ptrs[chain] as usize + i - 1] };
                        &codes[lcs(source).min(lcs_cap) as usize]
                    },
                };
                let len = match code.read(&mut bits).ok_or_else(|| truncated("suffix lengths"))? {
                    symbol if symbol == long_symbol => {
                        let long = bits.read_gamma().ok_or_else(|| truncated("suffix lengths"))? + max_direct as u64;
                        long.try_into().map_err(|_| truncated("suffix lengths"))?
                    },
                    symbol if symbol < long_symbol => symbol as u32,
                    _ => return Err(truncated("suffix lengths")),
                };
                suffix_lens.push(len);
            }
        }

        let mut targets = Vec::with_capacity(n_chains);
        let mut bits = BitReader::new(&stored.targets);
        let mut next_new = 2_u32;
        for _ in 0..n_chains {
            if bits.read_bits(1).ok_or_else(|| truncated("targets"))? == 1 {
                targets.push(next_new);
                next_new += 1;
            } else {
                targets.push(bits.read_truncated_binary(next_new).ok_or_else(|| truncated("targets"))?);
            }
        }
        if next_new as usize != n_junctions {
            return Err(truncated("targets"))
        }

        Ok(ChainCsr {
            junction_weights,
            row_ptrs,
            targets,
            chain_ptrs,
            inner_weights,
            suffix_lens,
        })
    }

    /// Number of linear chains.
    pub fn n_chains(&self) -> usize {
        self.targets.len()
    }

    /// Range of chains leaving junction `junction`.
    pub fn row(
        &self,
        junction: usize,
    ) -> (usize, usize) {
        (self.row_ptrs[junction] as usize, self.row_ptrs[junction + 1] as usize)
    }

    /// Colex ranks of the inner nodes of `chain`.
    pub fn inner_nodes(
        &self,
        chain: usize,
    ) -> &[u32] {
        &self.inner_weights[self.chain_ptrs[chain] as usize..self.chain_ptrs[chain + 1] as usize]
    }

    /// Suffix lengths of the edges in `chain`, one more than its inner nodes.
    pub fn chain_suffix_lens(
        &self,
        chain: usize,
    ) -> &[u32] {
        let start = self.chain_ptrs[chain] as usize + chain;
        let end = self.chain_ptrs[chain + 1] as usize + chain + 1;
        &self.suffix_lens[start..end]
    }

    /// Expands the chains back into an uncollapsed graph.
    ///
    /// Junctions keep their indices and inner nodes are numbered after them.
    /// Edges out of a junction keep the order of the chains, so the choice
    /// stream stays valid for the expanded graph.
    pub fn to_csr(
        &self,
    ) -> Csr<u32> {
        let n_junctions = self.junction_weights.len() as u32;
        let mut node_weights = self.junction_weights.clone();
        node_weights.extend(self.inner_weights.iter().copied());

        // (source, target, suffix length); junction rows come first in chain order.
        let mut edges: Vec<(u32, u32, u32)> = Vec::with_capacity(self.suffix_lens.len());
        for junction in 0..n_junctions as usize {
            let (start, end) = self.row(junction);
            for chain in start..end {
                let suffix_lens = self.chain_suffix_lens(chain);
                let mut prev = junction as u32;
                for (i, suffix_len) in suffix_lens.iter().enumerate().take(self.inner_nodes(chain).len()) {
                    let node = n_junctions + self.chain_ptrs[chain] + i as u32;
                    edges.push((prev, node, *suffix_len));
                    prev = node;
                }
                edges.push((prev, self.targets[chain], *suffix_lens.last().unwrap()));
            }
        }
        edges.sort_by_key(|edge| edge.0);

        let mut row_ptrs = vec![0_u32; node_weights.len() + 1];
        for edge in &edges {
            row_ptrs[edge.0 as usize + 1] += 1;
        }
        for i in 0..node_weights.len() {
            row_ptrs[i + 1] += row_ptrs[i];
        }

        Csr {
            node_weights,
            row_ptrs,
            column_indices: edges.iter().map(|edge| edge.1).collect(),
            suffix_lens: edges.iter().map(|edge| edge.2).collect(),
        }
    }
}

impl Walkable for ChainCsr {
    /// Walks the path of `color` from the source to the sink.
    ///
    /// Same as the walk of the uncollapsed graph; the bound from [Walkable]
    /// counts uncollapsed edges.
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError> {
        let mut reader = choices
            .reader(color)
            .ok_or_else(|| GraphDecodeError(format!("no choice stream for colour {}", color)))?;
        let mut junction = SOURCE_NODE as usize;
        let mut sequence_length = 0_u32;
        let mut dictionary: Vec<(u32, u32)> = Vec::new();
        let mut n_edges = 0_usize;

        loop {
            let chain = take_edge(self.row(junction), &mut reader, color, junction)?;
            let inner = self.inner_nodes(chain);
            let suffix_lens = self.chain_suffix_lens(chain);

            if junction == SOURCE_NODE as usize {
                sequence_length = suffix_lens[0];
            }
            n_edges += suffix_lens.len();
            let max_edges = (sequence_length as usize).saturating_add(1);
            if n_edges > max_edges {
                return Err(GraphDecodeError(format!("colour {} did not reach the sink in {} edges", color, max_edges)))
            }

            if junction != SOURCE_NODE as usize {
                dictionary.push((self.junction_weights[junction], suffix_lens[0]));
            }
            dictionary.extend(inner.iter().copied().zip(suffix_lens[1..].iter().copied()));

            junction = self.targets[chain] as usize;
            if junction == SINK_NODE as usize {
                if !reader.is_finished() {
                    return Err(GraphDecodeError(format!("unused choices left for colour {}", color)))
                }
                return Ok((sequence_length, dictionary))
            }
        }
    }
}

impl serde::Serialize for ChainCsr {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.to_stored().serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for ChainCsr {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let stored = StoredChains::deserialize(deserializer)?;
        ChainCsr::from_stored(stored).map_err(serde::de::Error::custom)
    }
}

