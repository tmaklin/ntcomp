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
//! The index: an SBWT of a _k_-mer set with the reverse complements and
//! select support, and its LCS array. Built, loaded and saved here; walked
//! through [IndexWalk].
//!
//! On disk, an index with prefix `p` is `p.sbwt` (the SBWT after the
//! length and name of its type, `SubsetMatrix`) and `p.lcs`.
//!
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

use sbwt::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::walk::{IndexWalk, WalkIndex};

type E = Box<dyn std::error::Error>;

/// The type name stored before the SBWT.
const SBWT_TYPE: &[u8] = b"SubsetMatrix";

/// Options for building an index.
#[derive(Clone, Debug, PartialEq)]
pub struct BuildOpts {
    /// _k_-mer size.
    pub k: usize,
    /// Threads.
    pub num_threads: usize,
    /// Length of the precalculated prefixes in the SBWT.
    pub prefix_precalc: usize,
    /// Memory for building on disk, in gigabytes.
    pub mem_gb: usize,
    /// Deduplicate the _k_-mer batches while building (saves memory).
    pub dedup_batches: bool,
    /// Build on disk in this directory instead of in memory.
    pub temp_dir: Option<String>,
}

impl Default for BuildOpts {
    fn default() -> Self {
        BuildOpts { k: 31, num_threads: 1, prefix_precalc: 8, mem_gb: 4, dedup_batches: false, temp_dir: None }
    }
}

/// An SBWT index with its LCS array, see the module docs.
pub struct Index {
    sbwt: SbwtIndex<SubsetMatrix>,
    lcs: LcsArray,
}

impl Index {
    /// The index of the _k_-mers of `sequences` and their reverse
    /// complements.
    pub fn build(
        sequences: &[Vec<u8>],
        opts: &BuildOpts,
    ) -> Self {
        assert!(!sequences.is_empty(), "no sequences to index");
        let (sbwt, lcs) = match &opts.temp_dir {
            Some(temp_dir) => BitPackedKmerSortingDisk::new_from_vecs(sequences, opts.k)
                .n_threads(opts.num_threads)
                .add_rev_comp(true)
                .build_lcs(true)
                .build_select_support(true)
                .precalc_length(opts.prefix_precalc)
                .mem_gb(opts.mem_gb)
                .dedup_batches(opts.dedup_batches)
                .temp_dir(PathBuf::from(OsString::from(temp_dir)).as_path())
                .run(),
            None => BitPackedKmerSortingMem::new_from_vecs(sequences, opts.k)
                .k(opts.k)
                .n_threads(opts.num_threads)
                .add_rev_comp(true)
                .build_lcs(true)
                .build_select_support(true)
                .precalc_length(opts.prefix_precalc)
                .dedup_batches(opts.dedup_batches)
                .run(),
        };
        Index { sbwt, lcs: lcs.expect("the LCS array is built") }
    }

    /// Loads the index `<prefix>.sbwt` and `<prefix>.lcs`.
    pub fn load(
        prefix: &str,
    ) -> Result<Self, E> {
        let open = |path: &str| File::open(path).map(BufReader::new).map_err(|e| format!("{}: {}", path, e));
        let sbwt_path = format!("{}.sbwt", prefix);
        let mut input = open(&sbwt_path)?;
        let invalid = |e: &dyn std::fmt::Display| format!("{}: not an SBWT ({})", sbwt_path, e);
        let mut type_len = [0_u8; 8];
        input.read_exact(&mut type_len).map_err(|e| invalid(&e))?;
        let mut type_name = vec![0_u8; (u64::from_le_bytes(type_len) as usize).min(64)];
        input.read_exact(&mut type_name).map_err(|e| invalid(&e))?;
        if type_name != SBWT_TYPE {
            return Err(invalid(&format!("type {}", String::from_utf8_lossy(&type_name))).into())
        }
        let sbwt = SbwtIndex::<SubsetMatrix>::load(&mut input).map_err(|e| invalid(&e))?;

        let lcs_path = format!("{}.lcs", prefix);
        let lcs = LcsArray::load(&mut open(&lcs_path)?).map_err(|e| format!("{}: not an LCS array ({})", lcs_path, e))?;
        if lcs.len() != sbwt.n_sets() {
            return Err(format!("{}: the LCS array has {} values for the {} rows of {}", lcs_path, lcs.len(), sbwt.n_sets(), sbwt_path).into())
        }
        Ok(Index { sbwt, lcs })
    }

    /// Saves the index as `<prefix>.sbwt` and `<prefix>.lcs`.
    pub fn save(
        &self,
        prefix: &str,
    ) -> Result<(), E> {
        let create = |path: &str| File::create(path).map(BufWriter::new).map_err(|e| format!("{}: {}", path, e));
        let sbwt_path = format!("{}.sbwt", prefix);
        let mut out = create(&sbwt_path)?;
        let write = |out: &mut BufWriter<File>| -> std::io::Result<()> {
            out.write_all(&(SBWT_TYPE.len() as u64).to_le_bytes())?;
            out.write_all(SBWT_TYPE)?;
            self.sbwt.serialize(out)?;
            out.flush()
        };
        write(&mut out).map_err(|e| format!("{}: {}", sbwt_path, e))?;

        let lcs_path = format!("{}.lcs", prefix);
        let mut out = create(&lcs_path)?;
        self.lcs.serialize(&mut out).and_then(|_| out.flush()).map_err(|e| format!("{}: {}", lcs_path, e))?;
        Ok(())
    }

    /// The _k_-mer size.
    pub fn k(&self) -> usize {
        self.sbwt.k()
    }

    /// Number of _k_-mers (with the reverse complements).
    pub fn n_kmers(&self) -> usize {
        self.sbwt.n_kmers()
    }

    /// The SBWT.
    pub fn sbwt(&self) -> &SbwtIndex<SubsetMatrix> {
        &self.sbwt
    }

    /// The LCS array.
    pub fn lcs(&self) -> &LcsArray {
        &self.lcs
    }

    /// For the matching statistics of sequences against the index.
    pub fn streaming(&self) -> StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray> {
        StreamingIndex::new(&self.sbwt, &self.lcs)
    }

    /// The index as a [WalkIndex] (with its LCS array).
    pub fn walk_index(&self) -> WalkIndex<'_> {
        WalkIndex { sbwt: &self.sbwt, lcs: Some(&self.lcs) }
    }

    /// The shortest match length that a random match against the index
    /// reaches with probability at most `max_error_prob` (at most _k_).
    pub fn random_match_threshold(
        &self,
        max_error_prob: f64,
    ) -> usize {
        random_match_threshold(self.k(), self.n_kmers(), 4, max_error_prob)
    }
}

impl IndexWalk for Index {
    fn sbwt(&self) -> &SbwtIndex<SubsetMatrix> {
        &self.sbwt
    }

    fn lcs_at(&self, row: usize) -> usize {
        self.lcs.access(row)
    }
}

/// Log of the CDF of the longest of `n_kmers` random matches of
/// `alphabet_size` letters at `t` (Mäklin et al. 2024, the kbo crate).
fn log_rm_max_cdf(
    t: usize,
    alphabet_size: usize,
    n_kmers: usize,
) -> f64 {
    n_kmers as f64 * (- ((1.0_f64.ln() - (alphabet_size as f64).ln()).exp()).powi(t as i32 + 1)).ln_1p()
}

/// The shortest match length `t` < `k` whose probability as a random
/// match against `n_kmers` _k_-mers is at most `max_error_prob`, else `k`.
fn random_match_threshold(
    k: usize,
    n_kmers: usize,
    alphabet_size: usize,
    max_error_prob: f64,
) -> usize {
    assert!(k > 0 && n_kmers > 0 && alphabet_size > 0);
    assert!(max_error_prob > 0_f64 && max_error_prob <= 1_f64);
    for i in 1..k {
        if log_rm_max_cdf(i, alphabet_size, n_kmers) > (-max_error_prob).ln_1p() {
            return i
        }
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::related_index;

    // The values of kbo's tests.
    #[test]
    fn random_match_cdf() {
        let expected = [-1306319.1078024083,-318761.2492719044,-79220.9269610741,-19776.1823255263,-4942.2344281681,-1235.4454790664,-308.8543003470,-77.2131332649,-19.3032557026,-4.8258121998,-1.2064529421,-0.3016132288,-0.0754033068,-0.0188508267,-0.0047127067,-0.0011781767,-0.0002945442,-0.0000736360,-0.0000184090,-0.0000046023,-0.0000011506,-0.0000002876,-0.0000000719,-0.0000000180,-0.0000000045,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000,0.0000000000];
        for t in 1..32 {
            assert!((log_rm_max_cdf(t, 4, 20240921) - expected[t - 1]).abs() < 1e-8, "t {}", t);
        }
    }

    #[test]
    fn random_match_thresholds() {
        for (i, expected) in [15, 18, 22, 25, 28].into_iter().enumerate() {
            assert_eq!(random_match_threshold(31, 20240921, 4, 0.01_f64.powf(i as f64 + 1.0)), expected);
        }
    }

    #[test]
    fn indexes_roundtrip_and_bad_files_are_errors() {
        let (_, index) = related_index(15);
        let dir = std::env::temp_dir().join(format!("ntcomp_index_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("index").to_string_lossy().to_string();
        index.save(&prefix).unwrap();
        let loaded = Index::load(&prefix).unwrap();
        assert_eq!((loaded.k(), loaded.n_kmers()), (15, index.n_kmers()));
        assert_eq!(loaded.lcs(), index.lcs());
        assert!(*loaded.sbwt() == *index.sbwt());

        // A missing or truncated LCS array, and one of another index.
        let other = dir.join("other").to_string_lossy().to_string();
        related_index(31).1.save(&other).unwrap();
        std::fs::copy(format!("{}.lcs", other), format!("{}.lcs", prefix)).unwrap();
        assert!(Index::load(&prefix).is_err());
        let bytes = std::fs::read(format!("{}.lcs", other)).unwrap();
        std::fs::write(format!("{}.lcs", prefix), &bytes[..bytes.len() / 2]).unwrap();
        assert!(Index::load(&prefix).is_err());
        std::fs::remove_file(format!("{}.lcs", prefix)).unwrap();
        assert!(Index::load(&prefix).is_err());
        // An empty SBWT file.
        std::fs::write(format!("{}.sbwt", prefix), b"").unwrap();
        assert!(Index::load(&prefix).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
