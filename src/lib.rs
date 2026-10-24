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
//! ntcomp: compressing assemblies against an SBWT index.
//!
//! - [archive]: compressing collections of assemblies into an archive and
//!   decompressing any of their contigs ([ArchiveWriter], [ArchiveReader]);
//! - [nlz]: the NLZ dictionaries of sequences against the index, and the
//!   plain NLZ file format;
//! - [graph]: NLZ-OG, the dictionaries of many sequences as coloured paths in
//!   a shared graph (the archive format);
//! - [tigs]: the index tigs, the _k_-mer set of the index;
//! - [index]: the SBWT index with its LCS array ([Index]); [walk]: walking
//!   its de Bruijn graph; [edits]: SNPs and
//!   small indels towards the index; [strands]: choosing the strand of each
//!   contig; [codecs]: the codecs shared by the formats; [fastx]: reading
//!   and writing FASTA and FASTQ files.
//!
pub mod archive;
pub mod codecs;
pub mod edits;
pub mod fastx;
pub mod graph;
pub mod index;
pub mod nlz;
pub mod strands;
pub mod walk;
pub mod tigs;
#[cfg(feature = "graph-view")]
pub mod graphview;

pub use archive::{ArchiveReader, ArchiveWriter};
pub use index::Index;

/// Test data shared by the modules' tests.
#[cfg(test)]
pub(crate) mod testing {
    use crate::index::{BuildOpts, Index};

    /// A random sequence of `len` bases from `seed`.
    pub(crate) fn random_sequence(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len).map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize]
        }).collect()
    }

    /// A random sequence of `n` bases from `seed`, by a xorshift generator
    /// (the tig tests' data).
    pub(crate) fn xorshift_sequence(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
        (0..n).map(|_| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b"ACGT"[(x >> 33) as usize % 4] }).collect()
    }

    /// Four related sequences (one with substitutions every 97, 211 and 401
    /// bases), so the de Bruijn graph branches, and their index for `k`.
    pub(crate) fn related_index(k: usize) -> (Vec<Vec<u8>>, Index) {
        let base = random_sequence(3000, 1);
        let mut seqs = vec![base.clone()];
        for (seed, step) in [(2, 97), (3, 211), (4, 401)] {
            let mut variant = base.clone();
            let flips = random_sequence(variant.len() / step + 1, seed);
            for (n, i) in (step / 2..variant.len()).step_by(step).enumerate() {
                variant[i] = flips[n];
            }
            seqs.push(variant);
        }
        let index = Index::build(&seqs, &BuildOpts { k, ..Default::default() });
        (seqs, index)
    }

}
