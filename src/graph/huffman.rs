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
//! Canonical Huffman codes.
//!
//! A code is stored as the code length of every symbol (0 for unused symbols);
//! the codes themselves are assigned canonically, in order of (code length,
//! symbol).
//!
use std::cmp::Reverse;
use std::collections::BinaryHeap;

use super::bits::{BitReader, BitWriter};

/// Longest code that fits in [BitWriter::write_bits].
const MAX_CODE_LEN: usize = 64;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HuffmanCode {
    /// Code length of every symbol, 0 if the symbol is not used.
    lengths: Vec<u8>,
    /// (code, code length) of every symbol.
    codes: Vec<(u64, u32)>,
    /// Number of codes of each length.
    counts: Vec<u64>,
    /// Used symbols in canonical order.
    symbols: Vec<u32>,
}

impl HuffmanCode {
    /// Builds a Huffman code for symbols `0..counts.len()` occurring `counts` times.
    pub fn from_counts(
        counts: &[u64],
    ) -> Self {
        let n = counts.len();
        let mut lengths = vec![0_u8; n];
        let used: Vec<usize> = (0..n).filter(|&s| counts[s] > 0).collect();
        if used.len() == 1 {
            lengths[used[0]] = 1;
        } else if used.len() > 1 {
            // Leaves are 0..n, internal nodes are added after them.
            let mut parent = vec![usize::MAX; n];
            let mut heap: BinaryHeap<Reverse<(u64, usize)>> = used.iter().map(|&s| Reverse((counts[s], s))).collect();
            while heap.len() > 1 {
                let Reverse((a, x)) = heap.pop().unwrap();
                let Reverse((b, y)) = heap.pop().unwrap();
                let node = parent.len();
                parent.push(usize::MAX);
                parent[x] = node;
                parent[y] = node;
                heap.push(Reverse((a + b, node)));
            }
            for &s in &used {
                let (mut depth, mut node) = (0, s);
                while parent[node] != usize::MAX {
                    node = parent[node];
                    depth += 1;
                }
                assert!(depth <= MAX_CODE_LEN, "Huffman code longer than {} bits", MAX_CODE_LEN);
                lengths[s] = depth as u8;
            }
        }
        Self::from_lengths(lengths).unwrap()
    }

    /// Rebuilds a code from the code lengths returned by [HuffmanCode::lengths].
    pub fn from_lengths(
        lengths: Vec<u8>,
    ) -> Result<Self, String> {
        let max_len = lengths.iter().copied().max().unwrap_or(0) as usize;
        if max_len > MAX_CODE_LEN {
            return Err(format!("Huffman code longer than {} bits", MAX_CODE_LEN))
        }

        let mut counts = vec![0_u64; max_len + 1];
        for &len in lengths.iter().filter(|&&len| len > 0) {
            counts[len as usize] += 1;
        }
        let kraft: u128 = (1..=max_len).map(|len| (counts[len] as u128) << (MAX_CODE_LEN - len)).sum();
        if kraft > 1_u128 << MAX_CODE_LEN {
            return Err("Huffman code lengths are not a prefix code".to_string())
        }

        let mut symbols: Vec<u32> = (0..lengths.len() as u32).filter(|&s| lengths[s as usize] > 0).collect();
        symbols.sort_by_key(|&s| (lengths[s as usize], s));

        let mut codes = vec![(0_u64, 0_u32); lengths.len()];
        let (mut code, mut prev_len) = (0_u64, 0_u32);
        for (i, &s) in symbols.iter().enumerate() {
            let len = lengths[s as usize] as u32;
            if i > 0 {
                code += 1;
            }
            code <<= len - prev_len;
            prev_len = len;
            codes[s as usize] = (code, len);
        }

        Ok(HuffmanCode { lengths, codes, counts, symbols })
    }

    /// Code length of every symbol, enough to rebuild the code.
    pub fn lengths(&self) -> &[u8] {
        &self.lengths
    }

    pub fn write(
        &self,
        symbol: usize,
        writer: &mut BitWriter,
    ) {
        let (code, len) = self.codes[symbol];
        assert!(len > 0, "symbol {} is not in the code", symbol);
        writer.write_bits(code, len);
    }

    /// Reads one symbol, None if the input ends or holds no valid code.
    pub fn read(
        &self,
        reader: &mut BitReader,
    ) -> Option<usize> {
        // Canonical decoding: codes of length `len` are first..first + counts[len].
        let (mut code, mut first, mut index) = (0_u64, 0_u64, 0_u64);
        for len in 1..self.counts.len() {
            code |= reader.read_bits(1)?;
            let count = self.counts[len];
            if code < first + count {
                return Some(self.symbols[(index + code - first) as usize] as usize)
            }
            index += count;
            first = (first + count) << 1;
            code <<= 1;
        }
        None
    }
}
