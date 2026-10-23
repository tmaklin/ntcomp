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
//! Choosing the strand of each contig.
//!
//! Assemblies of related isolates often have homologous contigs on opposite
//! strands. A contig and its reverse complement have different _k_-mers, so
//! the graph cannot share them. Storing every contig on the strand that
//! shares the most with the contigs before it lets the graph share them.
//!
//! The strands are chosen in input order from content-sampled _k_-mers: the
//! [SKETCH_K]-mers of a strand whose hash is 0 modulo the sampling width
//! (1 in [DEFAULT_SAMPLING] by default). Homologous contigs sample the same
//! _k_-mers wherever they are in the contig. A contig is reversed if more of
//! the sampled _k_-mers of its reverse complement than of its forward strand
//! are among the sampled _k_-mers of the strands chosen so far.
//!
use std::collections::HashSet;
use std::hash::{BuildHasherDefault, Hasher};

/// Length of the sampled _k_-mers.
pub const SKETCH_K: usize = 31;

/// Default sampling width: one in this many _k_-mers is sampled.
pub const DEFAULT_SAMPLING: u64 = 32;

/// Reverse complement of a nucleotide sequence; other characters are kept.
pub fn reverse_complement(
    seq: &[u8],
) -> Vec<u8> {
    seq.iter().rev().map(|c| match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => *other,
    }).collect()
}

/// splitmix64 finalizer.
fn mix(
    mut x: u64,
) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

/// Hashes of the [SKETCH_K]-mers of `seq` that are 0 modulo `sampling`. No
/// sampled _k_-mer contains a character other than A, C, G or T.
fn sketch(
    seq: &[u8],
    sampling: u64,
) -> Vec<u64> {
    let mask = (1_u64 << (2 * SKETCH_K)) - 1;
    let (mut kmer, mut len) = (0_u64, 0_usize);
    let mut hashes = Vec::new();
    for &c in seq {
        let bits = match c {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => { len = 0; continue },
        };
        kmer = ((kmer << 2) | bits) & mask;
        len += 1;
        if len >= SKETCH_K {
            let hash = mix(kmer);
            if hash % sampling == 0 {
                hashes.push(hash);
            }
        }
    }
    hashes
}

/// Hasher for values that are already hashes.
#[derive(Default)]
struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 << 8) | b as u64;
        }
    }

    fn write_u64(&mut self, x: u64) {
        self.0 = x;
    }
}

/// Chooses the strand of each contig in input order.
pub struct StrandChooser {
    /// One in `sampling` _k_-mers is sampled.
    sampling: u64,

    /// Sampled _k_-mers of the strands chosen so far.
    seen: HashSet<u64, BuildHasherDefault<IdentityHasher>>,
}

impl StrandChooser {
    /// Samples one in `sampling` _k_-mers; panics if `sampling` is 0.
    pub fn new(
        sampling: u64,
    ) -> Self {
        assert!(sampling > 0, "sampling width must be at least 1");
        StrandChooser { sampling, seen: HashSet::default() }
    }

    /// Returns true if `sequence` should be stored reverse complemented, and
    /// remembers the sampled _k_-mers of the chosen strand.
    pub fn choose(
        &mut self,
        sequence: &[u8],
    ) -> bool {
        let sketches = self.sketch(sequence);
        self.choose_sketched(sketches)
    }

    /// Sampled _k_-mers of both strands of `sequence`, for
    /// [StrandChooser::choose_sketched].
    ///
    /// Does not depend on the strands chosen so far, so the sketches of
    /// several contigs can be computed in parallel.
    pub fn sketch(
        &self,
        sequence: &[u8],
    ) -> StrandSketches {
        StrandSketches {
            forward: sketch(sequence, self.sampling),
            reverse: sketch(&reverse_complement(sequence), self.sampling),
        }
    }

    /// Same as [StrandChooser::choose] for a contig sketched with
    /// [StrandChooser::sketch]. Contigs must be chosen in input order.
    pub fn choose_sketched(
        &mut self,
        sketches: StrandSketches,
    ) -> bool {
        let StrandSketches { forward, reverse } = sketches;
        let shared = |hashes: &[u64]| hashes.iter().filter(|h| self.seen.contains(h)).count();
        let reversed = shared(&reverse) > shared(&forward);
        self.seen.extend(if reversed { reverse } else { forward });
        reversed
    }
}

/// Sampled _k_-mers of both strands of a contig, see [StrandChooser::sketch].
pub struct StrandSketches {
    forward: Vec<u64>,
    reverse: Vec<u64>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_sequence(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len).map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize]
        }).collect()
    }

    #[test]
    fn reverse_complement_roundtrip() {
        assert_eq!(reverse_complement(b"AACGTN"), b"NACGTT".to_vec());
        let seq = random_sequence(1000, 1);
        assert_eq!(reverse_complement(&reverse_complement(&seq)), seq);
    }

    #[test]
    fn homologous_contigs_get_the_same_strand() {
        let first = random_sequence(20000, 2);
        // A variant of the first contig on the opposite strand, and an unrelated contig.
        let mut variant = first.clone();
        for i in (500..variant.len()).step_by(2000) {
            variant[i] = if variant[i] == b'A' { b'C' } else { b'A' };
        }
        let flipped = reverse_complement(&variant);
        let unrelated = random_sequence(20000, 3);

        let mut chooser = StrandChooser::new(DEFAULT_SAMPLING);
        assert!(!chooser.choose(&first), "the first contig is kept as given");
        assert!(chooser.choose(&flipped), "the reverse complemented variant is flipped back");
        assert!(!chooser.choose(&variant), "the variant on the same strand is kept");
        assert!(!chooser.choose(&unrelated), "an unrelated contig is kept as given");

        // Sketching first and choosing later gives the same strands.
        let mut sketched = StrandChooser::new(DEFAULT_SAMPLING);
        let sketches: Vec<StrandSketches> = [&first, &flipped, &variant, &unrelated].iter().map(|seq| sketched.sketch(seq)).collect();
        let strands: Vec<bool> = sketches.into_iter().map(|s| sketched.choose_sketched(s)).collect();
        assert_eq!(strands, vec![false, true, false, false]);
    }

    #[test]
    fn short_and_ambiguous_contigs_are_kept() {
        let mut chooser = StrandChooser::new(DEFAULT_SAMPLING);
        assert!(!chooser.choose(b"ACGT"));
        assert!(!chooser.choose(&vec![b'N'; 1000]));
        assert!(!chooser.choose(b""));
    }
}
