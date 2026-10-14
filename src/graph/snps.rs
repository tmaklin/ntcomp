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
//! Editing SNPs towards the index before encoding.
//!
//! A contig that differs from the index by a SNP is encoded as several short
//! dictionary entries around the SNP. Replacing the contig's base with the
//! index's base before encoding gives fewer, longer entries; the original base
//! is stored separately and put back after decoding.
//!
//! The SNPs are called with [kbo::call] against the index. A SNP pays off if
//! few contigs share it: the graph already stores the nodes of a SNP that many
//! contigs share only once. SNPs are identified across contigs by the change
//! and [FLANK] bases on both sides.
//!
//! ## Stored SNPs
//!
//! Per colour, in contig order: the gap from the base after the previous SNP
//! (from the start of the contig for the first one) as a Huffman-coded bit
//! length followed by the lower bits, and the original base as its distance
//! from the index's base (1, 2 or 3 in A, C, G, T order; 2 is a transition),
//! also Huffman coded. The index's base is read from the decoded (edited)
//! contig. The streams of all colours are compressed together with LZMA.
//!
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use sbwt::LcsArray;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

use super::bits::{BitReader, BitWriter};
use super::huffman::HuffmanCode;
use super::lzma_wrapper::{self, LzmaData};

/// Default for the largest number of contigs that may share a SNP for it to
/// be edited.
pub const DEFAULT_MAX_SHARED: usize = 5;

/// Bases on both sides of a SNP that identify it across contigs.
pub const FLANK: usize = 32;

/// Number of gap symbols: the bit length of gap + 1, at most 64.
const N_GAP_SYMBOLS: usize = 65;

/// Number of base symbols: the distance from the index's base, 1..=3.
const N_BASE_SYMBOLS: usize = 4;

/// A SNP of a contig against the index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Snp {
    /// Position in the contig.
    pub pos: u32,

    /// The contig's base.
    pub base: u8,

    /// The index's base.
    pub index_base: u8,
}

/// A SNP and the hash that identifies it across contigs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CalledSnp {
    pub snp: Snp,
    pub key: u64,
}

fn base_code(
    base: u8,
) -> Option<u8> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Hash of the SNP with [FLANK] bases on both sides.
fn key(
    sequence: &[u8],
    snp: &Snp,
) -> u64 {
    let pos = snp.pos as usize;
    let start = pos.saturating_sub(FLANK);
    let end = (pos + 1 + FLANK).min(sequence.len());
    let mut hasher = DefaultHasher::new();
    sequence[start..pos].hash(&mut hasher);
    snp.base.hash(&mut hasher);
    snp.index_base.hash(&mut hasher);
    sequence[(pos + 1)..end].hash(&mut hasher);
    hasher.finish()
}

/// SNPs among the variants kbo calls for `variants`: sorted, without
/// overlapping calls, and with A, C, G or T on both sides.
fn snps_from_variants(
    sequence: &[u8],
    mut variants: Vec<kbo::variant_calling::Variant>,
) -> Vec<Snp> {
    variants.sort_by_key(|v| v.query_pos);
    let mut snps: Vec<Snp> = Vec::new();
    let mut prev_end = 0;
    let mut prev_pos = None;
    for v in variants {
        // Calls that overlap an earlier one are dropped, whatever their type.
        let end = v.query_pos + v.query_chars.len();
        if (v.query_chars.is_empty() && v.ref_chars.is_empty())
            || end > sequence.len()
            || v.query_pos < prev_end
            || prev_pos == Some(v.query_pos)
            || !v.ref_chars.iter().all(|&c| base_code(c).is_some()) {
            continue;
        }
        prev_end = end;
        prev_pos = Some(v.query_pos);
        if v.query_chars.len() == 1 && v.ref_chars.len() == 1 {
            let (base, index_base) = (sequence[v.query_pos], v.ref_chars[0]);
            if base_code(base).is_some() && base != index_base {
                snps.push(Snp { pos: v.query_pos.try_into().expect("contig position does not fit in u32"), base, index_base });
            }
        }
    }
    snps
}

/// Calls the SNPs of `sequence` against the index.
pub fn call_snps(
    sequence: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> Vec<CalledSnp> {
    let SbwtIndexVariant::SubsetMatrix(index) = sbwt;
    let mut opts = kbo::CallOpts::default();
    opts.sbwt_build_opts.k = index.k();
    opts.sbwt_build_opts.num_threads = 1;
    let variants = kbo::call(sbwt, lcs, sequence, opts);
    snps_from_variants(sequence, variants).into_iter().map(|snp| CalledSnp { snp, key: key(sequence, &snp) }).collect()
}

/// `sequence` with the index's base at every SNP.
pub fn edit(
    sequence: &[u8],
    snps: &[Snp],
) -> Vec<u8> {
    let mut edited = sequence.to_vec();
    for snp in snps {
        edited[snp.pos as usize] = snp.index_base;
    }
    edited
}

/// Error from reading the stored SNPs.
#[derive(Debug)]
pub struct SnpDecodeError(pub String);

impl std::fmt::Display for SnpDecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SNP decoding failed: {}", self.0)
    }
}

impl std::error::Error for SnpDecodeError {}

/// Edited SNPs of one colour.
#[derive(serde::Serialize, serde::Deserialize, Debug, Default, Clone, PartialEq, Eq)]
pub struct SnpStream {
    /// Number of SNPs in the stream.
    pub n_snps: u32,

    /// Huffman-coded gaps and bases, most significant bit first.
    pub bits: Vec<u8>,
}

/// Edited SNPs of every colour, see the module docs.
#[derive(serde::Serialize, serde::Deserialize, Debug, Default, Clone, PartialEq, Eq)]
struct StoredSnps {
    /// Huffman code lengths of the gap symbols.
    gap_code: Vec<u8>,

    /// Huffman code lengths of the base symbols.
    base_code: Vec<u8>,

    /// One stream per colour.
    streams: Vec<SnpStream>,
}

/// Edited SNPs of every colour.
///
/// Serialized as LZMA-compressed [StoredSnps].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnpEdits {
    gap_code: HuffmanCode,
    base_code: HuffmanCode,
    streams: Vec<SnpStream>,
}

fn gap_symbol(
    gap: u64,
) -> usize {
    (64 - (gap + 1).leading_zeros()) as usize
}

/// Distance of `base` from `index_base`, 1..=3.
fn base_symbol(
    base: u8,
    index_base: u8,
) -> usize {
    ((base_code(base).unwrap() + 4 - base_code(index_base).unwrap()) % 4) as usize
}

impl SnpEdits {
    /// Codes the SNPs of every colour, each sorted by position.
    pub fn new(
        snps: &[Vec<Snp>],
    ) -> Self {
        let mut gap_counts = vec![0_u64; N_GAP_SYMBOLS];
        let mut base_counts = vec![0_u64; N_BASE_SYMBOLS];
        for colour in snps {
            let mut next = 0_u64;
            for snp in colour {
                gap_counts[gap_symbol(snp.pos as u64 - next)] += 1;
                base_counts[base_symbol(snp.base, snp.index_base)] += 1;
                next = snp.pos as u64 + 1;
            }
        }
        let gap_code = HuffmanCode::from_counts(&gap_counts);
        let base_code = HuffmanCode::from_counts(&base_counts);

        let streams = snps.iter().map(|colour| {
            let mut bits = BitWriter::new();
            let mut next = 0_u64;
            for snp in colour {
                assert!(snp.pos as u64 >= next, "SNPs are not sorted by position");
                let gap = snp.pos as u64 - next;
                let symbol = gap_symbol(gap);
                gap_code.write(symbol, &mut bits);
                bits.write_bits(gap + 1, symbol as u32 - 1);
                base_code.write(base_symbol(snp.base, snp.index_base), &mut bits);
                next = snp.pos as u64 + 1;
            }
            SnpStream { n_snps: colour.len().try_into().expect("too many SNPs"), bits: bits.finish() }
        }).collect();

        SnpEdits { gap_code, base_code, streams }
    }

    /// Number of colours.
    pub fn n_colors(&self) -> usize {
        self.streams.len()
    }

    /// Total number of edited SNPs.
    pub fn n_snps(&self) -> usize {
        self.streams.iter().map(|s| s.n_snps as usize).sum()
    }

    /// Puts the original bases of `color` back into its decoded contig.
    pub fn revert(
        &self,
        color: u32,
        sequence: &mut [u8],
    ) -> Result<(), SnpDecodeError> {
        let stream = self.streams.get(color as usize).ok_or_else(|| SnpDecodeError(format!("no SNP stream for colour {}", color)))?;
        let truncated = || SnpDecodeError(format!("SNP stream of colour {} is truncated", color));
        let mut bits = BitReader::new(&stream.bits);
        let mut next = 0_u64;
        for _ in 0..stream.n_snps {
            let symbol = self.gap_code.read(&mut bits).ok_or_else(truncated)?;
            if symbol == 0 {
                return Err(truncated())
            }
            let gap = ((1_u64 << (symbol - 1)) | bits.read_bits(symbol as u32 - 1).ok_or_else(truncated)?) - 1;
            let pos = next + gap;
            let base = self.base_code.read(&mut bits).ok_or_else(truncated)?;
            let index_base = sequence.get(pos as usize).and_then(|&b| base_code(b))
                .ok_or_else(|| SnpDecodeError(format!("SNP of colour {} at {} is not on an A, C, G or T of the contig", color, pos)))?;
            if base == 0 {
                return Err(SnpDecodeError(format!("SNP of colour {} at {} does not change the base", color, pos)))
            }
            sequence[pos as usize] = b"ACGT"[(index_base as usize + base) % 4];
            next = pos + 1;
        }
        Ok(())
    }

    fn to_stored(
        &self,
    ) -> Result<LzmaData, String> {
        let stored = StoredSnps {
            gap_code: self.gap_code.lengths().to_vec(),
            base_code: self.base_code.lengths().to_vec(),
            streams: self.streams.clone(),
        };
        let bytes = postcard::to_allocvec(&stored).map_err(|e| format!("SNPs: {}", e))?;
        lzma_wrapper::compress(&bytes).map_err(|e| format!("SNPs: {}", e))
    }

    fn from_stored(
        stored: LzmaData,
    ) -> Result<Self, String> {
        let bytes = lzma_wrapper::decompress(&stored).map_err(|e| format!("SNPs: {}", e))?;
        let stored: StoredSnps = postcard::from_bytes(&bytes).map_err(|e| format!("SNPs: {}", e))?;
        if stored.gap_code.len() != N_GAP_SYMBOLS || stored.base_code.len() != N_BASE_SYMBOLS {
            return Err("SNPs: invalid code tables".to_string())
        }
        Ok(SnpEdits {
            gap_code: HuffmanCode::from_lengths(stored.gap_code).map_err(|e| format!("SNP gaps: {}", e))?,
            base_code: HuffmanCode::from_lengths(stored.base_code).map_err(|e| format!("SNP bases: {}", e))?,
            streams: stored.streams,
        })
    }
}

impl serde::Serialize for SnpEdits {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.to_stored().map_err(serde::ser::Error::custom)?.serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for SnpEdits {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let stored = LzmaData::deserialize(deserializer)?;
        SnpEdits::from_stored(stored).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kbo::variant_calling::Variant;

    fn snp(pos: u32, base: u8, index_base: u8) -> Snp {
        Snp { pos, base, index_base }
    }

    #[test]
    fn edit_and_revert_roundtrip() {
        let original = b"ACGTACGTACGTTTGCA".to_vec();
        let colours = vec![
            vec![snp(0, b'A', b'G'), snp(1, b'C', b'T'), snp(7, b'T', b'A'), snp(16, b'A', b'C')],
            vec![],
            vec![snp(5, b'C', b'G')],
        ];
        let edits: SnpEdits = postcard::from_bytes(&postcard::to_allocvec(&SnpEdits::new(&colours)).unwrap()).unwrap();
        assert_eq!(edits.n_colors(), 3);
        assert_eq!(edits.n_snps(), 5);
        for (color, snps) in colours.iter().enumerate() {
            let mut sequence = edit(&original, snps);
            if !snps.is_empty() {
                assert_ne!(sequence, original);
            }
            edits.revert(color as u32, &mut sequence).unwrap();
            assert_eq!(sequence, original);
        }
        assert!(edits.revert(3, &mut original.clone()).is_err());
    }

    #[test]
    fn corrupt_streams_are_errors() {
        let edits = SnpEdits::new(&[vec![snp(3, b'A', b'C'), snp(1000, b'G', b'A')]]);
        // The contig is too short for the second SNP.
        assert!(edits.revert(0, &mut b"ACGTACGT".to_vec()).is_err());
        let mut truncated = edits.clone();
        truncated.streams[0].bits.truncate(1);
        assert!(truncated.revert(0, &mut vec![b'A'; 2000]).is_err());
    }

    fn random_sequence(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len).map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize]
        }).collect()
    }

    #[test]
    fn builder_edits_snps_and_decoding_puts_them_back() {
        use super::super::{GraphBuilder, SuffixLen, Walkable, decode_path};

        let reference = random_sequence(4000, 7);
        // kbo calls SNPs only if k is well above the length of random matches.
        let k = 63;
        let mut opts = kbo::BuildOpts::default();
        opts.k = k as usize;
        opts.add_revcomp = true;
        opts.build_select = true;
        let (sbwt, lcs) = kbo::build(&[reference.clone()], opts);

        // Two contigs with the same SNPs and one with its own, more than k apart.
        let with_snps = |positions: &[usize]| {
            let mut contig = reference.clone();
            for &i in positions {
                contig[i] = b"ACGT"[(base_code(contig[i]).unwrap() as usize + 1) % 4];
            }
            contig
        };
        let contigs = [with_snps(&[500, 1500, 2500]), with_snps(&[500, 1500, 2500]), with_snps(&[1000, 3000])];

        for max_shared in [1, 2] {
            let mut builder = GraphBuilder::new().collapse_chains(true).edit_snps(true).edit_snps_max_shared(max_shared);
            assert!(builder.add_contig(&contigs[0], &sbwt, &lcs).is_err(), "SNPs must be called first");
            for contig in &contigs {
                builder.call_snps(contig, &sbwt, &lcs);
            }
            for contig in &contigs {
                builder.add_contig(contig, &sbwt, &lcs).unwrap();
            }
            let edits = builder.snp_edits().unwrap();
            // The shared SNPs are only edited if two contigs may share a SNP.
            assert_eq!(edits.n_snps(), if max_shared == 1 { 2 } else { 8 });
            let (graph, choices) = builder.finish();
            for (color, contig) in contigs.iter().enumerate() {
                let (len, dictionary) = graph.walk_path(&choices, color as u32).unwrap();
                let dictionary = dictionary.into_iter().map(|(colex, len)| (colex, if len == 0 {
                    SuffixLen::Implicit
                } else if len > k {
                    SuffixLen::Explicit { suffix_len: len }
                } else {
                    SuffixLen::Plain { suffix_len: len }
                })).collect();
                let mut decoded = decode_path(dictionary, len as usize, &sbwt, Some(&lcs)).unwrap();
                edits.revert(color as u32, &mut decoded).unwrap();
                assert_eq!(&decoded, contig, "colour {}", color);
            }
        }
    }

    #[test]
    fn only_non_overlapping_snps_are_kept() {
        let sequence = b"AAAAACCCCCGGGGGTTTTT".to_vec();
        let variant = |pos: usize, query: &[u8], reference: &[u8]| Variant { query_pos: pos, query_chars: query.to_vec(), ref_chars: reference.to_vec() };
        let variants = vec![
            variant(12, b"G", b"T"),   // SNP
            variant(2, b"AAA", b""),   // deletion from the index, covers 2..5
            variant(3, b"A", b"C"),    // overlaps the deletion
            variant(12, b"G", b"A"),   // second call at the same position
            variant(16, b"T", b"N"),   // not a base in the index
            variant(7, b"C", b"GG"),   // not a SNP
            variant(9, b"C", b"A"),    // SNP
        ];
        assert_eq!(snps_from_variants(&sequence, variants), vec![snp(9, b'C', b'A'), snp(12, b'G', b'T')]);
    }
}
