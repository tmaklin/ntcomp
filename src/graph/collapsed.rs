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
//! ## Stored suffix lengths
//!
//! The suffix lengths on the edges out of the dummy source node are the total
//! contig lengths and are stored as they are.
//!
//! Every other suffix length is a symbol of one canonical Huffman code, in
//! chain order:
//! - 0: the length is implicit, see [SuffixLen::Implicit](super::SuffixLen),
//! - 1..=255: the length itself,
//! - 256: a longer length, followed by the length - 255 (Elias gamma code).
//!
//! TODO Do not hardcode 255 as the maximum length -- this value should be taken
//! from the SBWT as it is equal to the _k_-mer length used.
//!
use super::{
    Csr,
    Choices,
    GraphDecodeError,
    SINK_NODE,
    SOURCE_NODE,
    take_edge,
    Walkable,
};
use super::bits::{BitReader, BitWriter};
use super::huffman::HuffmanCode;

/// Longest suffix length that is its own Huffman symbol.
///
/// This is equal to the value of `k` from the SBWT that was used to generate
/// the graph.
///
/// TODO Should be stored somewhere instead of being hardcoded here.
///
const MAX_DIRECT_SUFFIX_LEN: u32 = 255;
/// Huffman symbol for longer suffix lengths, followed by the length in Elias gamma.
///
/// TODO Infer from MAX_DIRECT_SUFFIX_LEN once it's stored inside this graph.
///
const LONG_SUFFIX_SYMBOL: usize = MAX_DIRECT_SUFFIX_LEN as usize + 1;

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
    /// Colex ranks of the junction nodes.
    pub junction_weights: Vec<u32>,

    /// Colex ranks of the linear chains, concatenated over all chains.
    pub inner_weights: Vec<u32>,

    /// Suffix lengths of the edges out of the dummy source node (the contig lengths).
    pub sequence_lens: Vec<u32>,

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
        let mut sequence_lens = Vec::with_capacity(n_source_chains);
        let mut symbols: Vec<(usize, u32)> = Vec::with_capacity(self.suffix_lens.len());
        for chain in 0..self.n_chains() {
            for (i, &len) in self.chain_suffix_lens(chain).iter().enumerate() {
                if chain < n_source_chains && i == 0 {
                    sequence_lens.push(len);
                } else if len > MAX_DIRECT_SUFFIX_LEN {
                    symbols.push((LONG_SUFFIX_SYMBOL, len));
                } else {
                    symbols.push((len as usize, len));
                }
            }
        }
        let mut counts = vec![0_u64; LONG_SUFFIX_SYMBOL + 1];
        for &(symbol, _) in &symbols {
            counts[symbol] += 1;
        }
        let code = HuffmanCode::from_counts(&counts);
        let mut suffix_lens = BitWriter::new();
        for &(symbol, len) in &symbols {
            code.write(symbol, &mut suffix_lens);
            if symbol == LONG_SUFFIX_SYMBOL {
                suffix_lens.write_gamma((len - MAX_DIRECT_SUFFIX_LEN) as u64);
            }
        }

        StoredChains {
            junction_weights: self.junction_weights.clone(),
            inner_weights: self.inner_weights.clone(),
            sequence_lens,
            suffix_code: code.lengths().to_vec(),
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
        let truncated = |what: &str| GraphDecodeError(format!("stored chains: {} truncated or invalid", what));
        let n_junctions = stored.junction_weights.len();
        if n_junctions < 2 {
            return Err(GraphDecodeError("stored chains: missing source or sink".to_string()))
        }

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
        if *chain_ptrs.last().unwrap() as usize != stored.inner_weights.len() {
            return Err(truncated("chain lengths"))
        }

        let n_source_chains = row_ptrs[1] as usize;
        if stored.sequence_lens.len() != n_source_chains {
            return Err(truncated("contig lengths"))
        }
        let code = HuffmanCode::from_lengths(stored.suffix_code).map_err(|e| GraphDecodeError(format!("stored chains: {}", e)))?;
        let mut bits = BitReader::new(&stored.suffix_lens);
        let mut suffix_lens = Vec::with_capacity(stored.inner_weights.len() + n_chains);
        for chain in 0..n_chains {
            let n_edges = (chain_ptrs[chain + 1] - chain_ptrs[chain]) as usize + 1;
            for i in 0..n_edges {
                if chain < n_source_chains && i == 0 {
                    suffix_lens.push(stored.sequence_lens[chain]);
                    continue;
                }
                let len = match code.read(&mut bits).ok_or_else(|| truncated("suffix lengths"))? {
                    LONG_SUFFIX_SYMBOL => {
                        let long = bits.read_gamma().ok_or_else(|| truncated("suffix lengths"))? + MAX_DIRECT_SUFFIX_LEN as u64;
                        long.try_into().map_err(|_| truncated("suffix lengths"))?
                    },
                    symbol if symbol < LONG_SUFFIX_SYMBOL => symbol as u32,
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
            junction_weights: stored.junction_weights,
            row_ptrs,
            targets,
            chain_ptrs,
            inner_weights: stored.inner_weights,
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
    /// Same as [super::walk_path] for the uncollapsed graph; `max_steps` bounds
    /// the number of uncollapsed edges.
    fn walk_path(
        &self,
        choices: &Choices,
        color: u32,
        max_steps: usize,
    ) -> Result<(u32, Vec<(u32, u32)>), GraphDecodeError> {
        let mut reader = choices
            .reader(color)
            .ok_or_else(|| GraphDecodeError(format!("no choice stream for colour {}", color)))?;
        let mut junction = SOURCE_NODE as usize;
        let mut sequence_length = 0_u32;
        let mut dictionary: Vec<(u32, u32)> = Vec::new();
        let mut steps = 0;

        loop {
            let chain = take_edge(self.row(junction), &mut reader, color, junction)?;
            let inner = self.inner_nodes(chain);
            let suffix_lens = self.chain_suffix_lens(chain);

            steps += suffix_lens.len();
            if steps > max_steps {
                return Err(GraphDecodeError(format!("colour {} did not reach the sink in {} steps", color, max_steps)))
            }

            if junction == SOURCE_NODE as usize {
                sequence_length = suffix_lens[0];
            } else {
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

