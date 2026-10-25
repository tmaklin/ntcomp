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
//! Stored out-degrees and chain lengths of [ChainCsr](super::ChainCsr).
//!
//! One arithmetic-coded stream ([crate::codecs::arith]) holds the out-degree of
//! every junction except the sink, in junction order, then the number of
//! inner nodes of every chain, in chain order, then the contig lengths (the
//! suffix lengths of the edges out of the dummy source, in chain order). Each
//! is an adaptive Elias-gamma code ([GammaModel](crate::codecs::arith::GammaModel): of the degree, of the number
//! of inner nodes + 1, of the contig length + 1).
//! The chain lengths have one model per context: the out-degree of the
//! chain's source junction (capped at 8) and the length of the previous chain
//! of the same junction (in five classes, or none): the chains of a junction
//! are alternative paths through the same region.
//!
use crate::codecs::arith::{Coder, GammaModel};
use crate::codecs::prealloc;
use crate::graph::{SINK_NODE, SOURCE_NODE};

/// Context of a chain length: the source junction's out-degree (capped at 8)
/// and the number of inner nodes of the previous chain of the junction.
fn chain_len_context(
    degree: u64,
    previous: Option<u64>,
) -> usize {
    let previous = match previous { None => 5, Some(0) => 0, Some(1) => 1, Some(2..=3) => 2, Some(4..=15) => 3, Some(_) => 4 };
    degree.min(8) as usize * 6 + previous
}

/// Out-degrees of the junctions (0 for the sink), the number of inner nodes of
/// every chain and the contig lengths, see [code_topology].
pub(super) type Topology = (Vec<u64>, Vec<u64>, Vec<u64>);

/// Codes the out-degrees of the junctions, the number of inner nodes of the
/// chains and the contig lengths (`shape`), or decodes them for `n_junctions`
/// junctions (`shape` None), see the module docs. None if decoded values are
/// out of range.
pub(super) fn code_topology<C: Coder>(
    c: &mut C,
    n_junctions: usize,
    shape: Option<(&[u64], &[u64], &[u64])>,
) -> Option<Topology> {
    let mut degree_model = GammaModel::new();
    let mut degrees = Vec::with_capacity(prealloc(n_junctions));
    let mut n_chains = 0_u64;
    for junction in 0..n_junctions {
        let degree = if junction == SINK_NODE as usize { 0 } else { degree_model.code(c, shape.map_or(1, |s| s.0[junction])) };
        degrees.push(degree);
        n_chains += degree;
        if n_chains > u32::MAX as u64 || c.overrun() { return None }
    }
    let mut models: Vec<GammaModel> = (0..9 * 6).map(|_| GammaModel::new()).collect();
    let mut inner = Vec::with_capacity(prealloc(n_chains as usize));
    for &degree in &degrees {
        let mut previous = None;
        for _ in 0..degree {
            let w = models[chain_len_context(degree, previous)].code(c, shape.map_or(1, |s| s.1[inner.len()] + 1));
            inner.push(w - 1);
            previous = Some(w - 1);
            if c.overrun() { return None }
        }
    }
    let mut length_model = GammaModel::new();
    let mut contig_lens = Vec::with_capacity(prealloc(degrees[SOURCE_NODE as usize] as usize));
    for i in 0..degrees[SOURCE_NODE as usize] as usize {
        let w = length_model.code(c, shape.map_or(1, |s| s.2[i] + 1));
        if w - 1 > u32::MAX as u64 || c.overrun() { return None }
        contig_lens.push(w - 1);
    }
    Some((degrees, inner, contig_lens))
}
