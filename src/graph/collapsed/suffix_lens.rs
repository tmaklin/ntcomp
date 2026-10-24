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
//! Stored suffix lengths of [ChainCsr](super::ChainCsr).
//!
//! The suffix lengths on the edges out of the dummy source node are the total
//! contig lengths and are stored as they are.
//!
//! Every other suffix length is coded in chain order with the binary
//! arithmetic coder of [crate::codecs::arith], with adaptive models per
//! context:
//! - a flag: the length is implicit (0, see
//!   [SuffixLen::Implicit](crate::nlz::SuffixLen::Implicit));
//! - a flag: the length is above D;
//! - a length up to D: length − 1 through an adaptive binary tree;
//! - a length above D: length − D as its bit length (a 6-bit adaptive tree
//!   shared by all contexts) and the bits below its top one, the first four
//!   adaptive given the bit length and the bits above them, the rest raw,
//!
//! where D = max(_k_, [MIN_DIRECT_SUFFIX_LEN]) (see [max_direct_suffix_len]),
//! with _k_ taken from the SBWT when the graph is stored or read with the
//! index, so every short entry (< _k_) is coded directly. Without the index D
//! is [MIN_DIRECT_SUFFIX_LEN], the same as for any _k_ ≤ 255. D is stored, so
//! reading with another _k_ fails.
//!
//! Without the LCS array there is one context for all edges. With it, the
//! context of an edge is min(LCS(r), D), where r is the colex rank of the
//! edge's source node: a short entry's rank is the start of the SBWT interval
//! of its suffix, so its length is always larger than LCS(r), often by one.
//! The decoder takes LCS(r) from the LCS array or computes it from the SBWT
//! with [crate::walk::lcs_from_sbwt].
//!
use crate::codecs::arith::{code_bit, code_tree, Coder, GammaModel, C0};

/// Smallest number of suffix lengths coded directly (not escaped): lengths
/// up to this are coded directly even when _k_ is smaller (left extensions
/// longer than _k_ code better directly than escaped).
pub const MIN_DIRECT_SUFFIX_LEN: u32 = 255;

/// Longest suffix length coded directly for _k_-mer size `k`
/// (None: unknown), see the module docs. Also the cap of the LCS contexts.
pub fn max_direct_suffix_len(
    k: Option<u32>,
) -> u32 {
    k.unwrap_or(0).max(MIN_DIRECT_SUFFIX_LEN)
}


/// Adaptive models of the suffix lengths (see the module docs) for
/// `n_contexts` contexts and direct lengths up to `max_direct`.
pub(super) struct SuffixModel {
    max_direct: u32,
    /// Bits of a direct length − 1.
    bits: u32,
    implicit: Vec<u32>,
    escaped: Vec<u32>,
    trees: Vec<u32>,
    long: GammaModel,
}

impl SuffixModel {
    pub(super) fn new(
        n_contexts: usize,
        max_direct: u32,
    ) -> Self {
        let bits = 32 - (max_direct - 1).leading_zeros();
        SuffixModel {
            max_direct, bits,
            implicit: vec![C0; n_contexts],
            escaped: vec![C0; n_contexts],
            trees: vec![C0; n_contexts << bits],
            long: GammaModel::new(),
        }
    }

    /// Codes (or decodes, `len` ignored) the suffix length of an edge in `context`.
    pub(super) fn code<C: Coder>(
        &mut self,
        c: &mut C,
        context: usize,
        len: u32,
    ) -> Option<u32> {
        if code_bit(c, &mut self.implicit[context], (len == 0) as u32) == 1 {
            return Some(0)
        }
        if code_bit(c, &mut self.escaped[context], (len > self.max_direct) as u32) == 0 {
            let tree = &mut self.trees[context << self.bits..(context + 1) << self.bits];
            let value = code_tree(c, tree, self.bits, len.wrapping_sub(1)) + 1;
            return (value <= self.max_direct).then_some(value)
        }
        let x = self.long.code(c, len.wrapping_sub(self.max_direct) as u64);
        u32::try_from(self.max_direct as u64).ok()?.checked_add(u32::try_from(x).ok()?)
    }
}
