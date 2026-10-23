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
use std::ops::Range;

use bincode::{Encode, Decode};
use bincode::encode_into_std_write;
use bincode::decode_from_slice;

use sbwt::LcsArray;
use sbwt::StreamingIndex;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

use walk::IndexWalk;

pub mod codecs;
pub mod encode;
pub mod decode;
pub mod edits;
pub mod graph;
pub mod strands;
pub mod walk;
pub mod tigs;
#[cfg(feature = "graph-view")]
pub mod graphview;

type E = Box<dyn std::error::Error>;

/// File header for encoded data
#[derive(Clone, Debug, Decode, Encode, PartialEq)]
pub struct FileHeader {
    /// Nlz header, consists of four ASCII bytes spelling "nrlz" and two bytes specifying the version.
    ///
    /// First four bytes can be used to check that a binary record is an .nlz record.
    /// Next two bytes can be used to check which version of nlz was used to generate this file.
    pub nlz_header: [u8; 6], // = [110, 114, 108, 122, ...];

    pub file_range_bytes: u32,

    pub hash_bytes: u32,

    /// Size of the section that marks the reverse complemented contigs.
    pub strand_bytes: u32,

    /// Size of the section with the edited SNPs.
    pub snp_bytes: u32,

    pub n_queries: u32,

    pub contig_name_bytes: u64,

    pub graph_bytes: u64,

    pub choice_bytes: u64,
}

pub const FILE_HEADER_BYTES: usize = 50;

#[derive(Encode, Decode)]
pub struct BlockHeader {
    pub block_size: u32,
    pub num_records: u32,
    pub num_u64: u32,
    pub encoded_size: u32,
    pub rice_param: u64,
    pub bitpacker_exponent: u8,
    pub placeholder1: u8,
    pub placeholder2: u32,
    pub placeholder3: u16,
    // TODO should store the used codec so others can be used
    // TODO store total number of bases?
}

pub fn encode_file_header(
    _ph1: u64,
    _ph2: u64,
    _ph3: u64,
    _ph4: u64
) -> Result<Vec<u8>, E> {
    let mut bytes: Vec<u8> = Vec::new();
    let header_placeholder = FileHeader{
        nlz_header: [0_u8; 6],
        strand_bytes: 0_u32,
        snp_bytes: 0_u32,
        n_queries: 0_u32,
        contig_name_bytes: 0_u64,
        graph_bytes: 0_u64,
        hash_bytes: 0_u32,
        file_range_bytes: 0_u32,
        choice_bytes: 0_u64,
    };
    let nbytes = encode_into_std_write(
        &header_placeholder,
        &mut bytes,
        bincode::config::standard().with_fixed_int_encoding(),
    )?;
    assert_eq!(nbytes, FILE_HEADER_BYTES);
    Ok(bytes)
}

pub fn decode_file_header(
    header_bytes: &[u8],
) -> Result<FileHeader, E> {
    Ok(decode_from_slice(header_bytes, bincode::config::standard().with_fixed_int_encoding())?.0)
}

pub fn encode_block_header(
    header: &BlockHeader,
) -> Result<Vec<u8>, E> {
    let mut bytes: Vec<u8> = Vec::new();
    let nbytes = encode_into_std_write(
        header,
        &mut bytes,
        bincode::config::standard().with_fixed_int_encoding(),
    )?;
    assert_eq!(nbytes, 32);
    Ok(bytes)
}

pub fn decode_block_header(
    header_bytes: &[u8],
) -> Result<BlockHeader, E> {
    Ok(decode_from_slice(header_bytes, bincode::config::standard().with_fixed_int_encoding())?.0)
}

/// Left extends `kmer_start` while the extension is unique in the SBWT, by at
/// most `max_extension_len` characters.
///
/// If `ref_nts` is given, also checks that the extension matches `ref_nts`.
///
/// Returns the extended k-mer and whether the extension stopped because the
/// SBWT has no unique left extension, i.e. where [left_extend_kmer] stops as
/// well. False if it stopped because of a mismatch with `ref_nts` or because
/// `max_extension_len` was reached.
///
/// Four searches per character: the reference the walks of
/// [IndexWalk::left_extend](walk::IndexWalk::left_extend) (one step per
/// character, from a known row) are checked against, and the left extension
/// of the plain NLZ format.
pub fn left_extend_kmer(
    kmer_start: &[u8],
    ref_nts: Option<&[u8]>,
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    max_extension_len: usize,
) -> (Vec<u8>, bool) {
    assert!(!kmer_start.is_empty());
    if let Some(ref_nts) = ref_nts {
        assert!(!ref_nts.is_empty());
    }

    let mut left_extension_len = 0;
    let mut kmer = kmer_start.to_vec();
    let mut stopped_in_sbwt = false;
    while left_extension_len < max_extension_len {
        let new_kmers: Vec<(Vec<u8>, Range<usize>)> = sbwt.alphabet().iter().filter_map(|c| {
            let new_kmer: Vec<u8> = [&[*c], &kmer[0..(kmer.len() - (left_extension_len + 1))]].concat();
            let res = sbwt.search(&new_kmer);
            if res.as_ref().is_some() {
                Some((new_kmer, res.unwrap()))
            } else {
                None
            }
        }).collect();
        if !new_kmers.is_empty() {
            let unique = new_kmers.len() == 1 && new_kmers[0].1.end - new_kmers[0].1.start == 1;

            let seq_matches = if let Some(ref_nts) = ref_nts {
                new_kmers[0].0[0] == ref_nts[ref_nts.len() - kmer.len() - 1]
            } else {
                true
            };

            if seq_matches && unique {
                kmer = [&[new_kmers[0].0[0]], kmer.as_slice()].concat();
            } else {
                stopped_in_sbwt = !unique;
                break;
            }
        } else {
            stopped_in_sbwt = true;
            break;
        }
        left_extension_len += 1;
    }
    (kmer, stopped_in_sbwt)
}

/// Encode a nucleotide sequence against an SBWT index.
///
/// Returns the dictionary and, for every dictionary entry, a boolean denoting
/// whether its suffix length is *implicit*: the entry is a left extended k-mer
/// whose extension stopped where [left_extend_kmer] stops as well, so the
/// decoder can recover the length by extending the k-mer as far as possible.
pub fn encode_sequence(
    nucleotides: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> Result<(Vec<(usize, Range<usize>)>, Vec<bool>), E> {
    let n = nucleotides.len();
    let mut implicit: Vec<bool> = Vec::new();
    let dictionary: Vec<(usize, Range<usize>)> = match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            let index = StreamingIndex::new(sbwt, lcs);
            let mut res = index.matching_statistics(nucleotides);
            let walk = walk::LcsIndex { sbwt, lcs };

            let mut i = n;
            let mut kept: Vec<usize> = Vec::new();
            while i > 0 {
                if res[i - 1].0 == k && i > k + 1 {
                    // Left extend as [left_extend_kmer] with the sequence as the
                    // reference, through the rows of the matching statistics: the
                    // k-mer ending one base to the left is the predecessor the
                    // sequence takes, and its group holds all the predecessors.
                    let max_extension_len = i - k - 1;
                    let mut row = res[i - 1].1.start;
                    let mut extension_len = 0;
                    let mut stopped_in_sbwt = false;
                    while extension_len < max_extension_len {
                        let left = i - 2 - extension_len;
                        if res[left].0 == k {
                            let predecessor = res[left].1.start;
                            if !walk.only_in_group(predecessor) {
                                stopped_in_sbwt = true;
                                break
                            }
                            row = predecessor;
                            extension_len += 1;
                        } else {
                            // The sequence's predecessor is not in the index: the
                            // extension stops in the SBWT unless the k-mer has exactly
                            // one (real) predecessor, a different one.
                            stopped_in_sbwt = !sbwt.inverse_lf_step(row).is_some_and(|predecessor| {
                                walk.only_in_group(predecessor) && walk.is_real_kmer(predecessor)
                            });
                            break
                        }
                    }
                    let phrase_len = k + extension_len;

                    // Store the whole left extension; the next entry ends where it starts.
                    res[i - 1].0 = phrase_len;
                    kept.push(i - 1);
                    implicit.push(stopped_in_sbwt);
                    i -= phrase_len - 1;
                } else {
                    kept.push(i - 1);
                    implicit.push(false);
                    if i > res[i - 1].0 {
                        i -= res[i - 1].0 - 1;
                    } else {
                        break;
                    }
                }

                if i > 0 {
                    i -= 1;
                } else {
                    break;
                }
            }
            kept.iter().map(|x| res[*x].clone()).collect()
        },
    };

    let dictionary_bases = dictionary.iter().map(|x| x.0).sum::<usize>();

    assert_eq!(dictionary_bases, nucleotides.len());
    assert_eq!(dictionary.len(), implicit.len());

    Ok((dictionary, implicit))
}

pub fn write_block_to<W: std::io::Write>(
    u64_encoding: &[u64],
    num_records: usize,
    sink: &mut W,
) -> Result<(), E> {
    let (data_1, data_2, data_3, data_4) = encode::split_encoded_dictionary(u64_encoding)?;

    assert_eq!(data_1.len(), data_2.len());
    assert_eq!(data_3.len(), u64_encoding.len());

    let block_1 = encode::compress_block(&data_1, num_records, encode::Codec::MinimalBinary)?;
    let block_2 = encode::compress_block(&data_2, num_records, encode::Codec::Rice)?;
    let block_3 = encode::compress_block(&data_3, num_records, encode::Codec::Rice)?;
    let block_4 = encode::compress_block(&data_4, num_records, encode::Codec::MinimalBinary)?;

    sink.write_all(&block_1)?;
    sink.write_all(&block_2)?;
    sink.write_all(&block_3)?;
    sink.write_all(&block_4)?;
    Ok(())
}

pub fn decode_sequence(
    encoding: &[u64],
    sbwt: &SbwtIndexVariant,
) -> Vec<Vec<u8>> {
    // TODO this needs to check if the sbwt has select support

    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut sequence: Vec<u8> = Vec::new();
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            let mut bases: usize = 0;
            encoding.iter().rev().for_each(|record| {
                let bytes: Vec<u8> = record.to_ne_bytes()[0..8].to_vec();
                let mut arr: [u8; 8] = [0; 8];
                arr[0..8].copy_from_slice(&bytes);
                let flags = bytes[7];
                let flag = u8::from_ne_bytes([flags]);
                let first: bool = (flag & 0b00000001) == 0b00000001;
                if flag & 0b00000010 == 0b00000000 {
                    let mut arr1: [u8; 4] = [0; 4];
                    let mut arr2: [u8; 4] = [0; 4];
                    let mut arr3: [u8; 1] = [0; 1];

                    arr1.copy_from_slice(&arr[0..4]);
                    arr2[0..3].copy_from_slice(&arr[4..7]);
                    arr3.copy_from_slice(&arr[7..8]);
                    let colex_rank = u32::from_ne_bytes(arr1);
                    let suffix_len = u32::from_ne_bytes(arr2);
                    bases += suffix_len as usize;

                    let kmer = if suffix_len > k as u32 {
                        let kmer = sbwt.access_kmer(colex_rank as usize);
                        let (new_kmer, _) = left_extend_kmer(&kmer, None, sbwt, (suffix_len - k as u32) as usize);
                        assert_eq!(new_kmer.len(), suffix_len as usize);
                        new_kmer
                    } else {
                        sbwt.access_kmer(colex_rank as usize)
                    };

                    sequence.extend(kmer[(kmer.len() - (suffix_len as usize))..kmer.len()].iter());
                } else {
                    let length: usize = ((flag & 0b11111100) >> 2) as usize;

                    let mut arr1: [u8; 8] = [0; 8];
                    arr1[0..7].copy_from_slice(&bytes[0..7]);

                    let bitnucs = u64::from_ne_bytes(arr1);
                    let mut kmer: Vec<u8> = Vec::new();
                    let _ = bitnuc::from_2bit(bitnucs, length, &mut kmer);
                    bases += kmer.len();
                    sequence.extend(kmer.iter());
                }

                if first {
                    bases = 0;
                    sequences.push(sequence.clone());
                    sequence.clear();
                }
            });
        },
    }

    sequences.into_iter().rev().collect()
}

pub fn decode_positions(
    encoding: &[u64],
    sbwt: &SbwtIndexVariant,
) -> Vec<(usize, usize, usize, usize, u8, Option<u32>)> {
    // TODO this needs to check if the sbwt has select support

    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut sequence: Vec<u8> = Vec::new();

    let mut positions: Vec<(usize, usize, usize, usize, u8, Option<u32>)> = Vec::new();
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            let mut bases: usize = 0;
            encoding.iter().rev().for_each(|record| {
                let start = bases;
                let encoding_type: u8;
                let colex: Option<u32>;

                let bytes: Vec<u8> = record.to_ne_bytes()[0..8].to_vec();
                let mut arr: [u8; 8] = [0; 8];
                arr[0..8].copy_from_slice(&bytes);
                let flags = bytes[7];
                let flag = u8::from_ne_bytes([flags]);
                let first: bool = (flag & 0b00000001) == 0b00000001;
                if flag & 0b00000010 == 0b00000000 {
                    let mut arr1: [u8; 4] = [0; 4];
                    let mut arr2: [u8; 4] = [0; 4];
                    let mut arr3: [u8; 1] = [0; 1];

                    arr1.copy_from_slice(&arr[0..4]);
                    arr2[0..3].copy_from_slice(&arr[4..7]);
                    arr3.copy_from_slice(&arr[7..8]);
                    let colex_rank = u32::from_ne_bytes(arr1);
                    let suffix_len = u32::from_ne_bytes(arr2);
                    bases += suffix_len as usize;

                    colex = Some(colex_rank);
                    let kmer = if suffix_len > k as u32 {
                        let kmer = sbwt.access_kmer(colex_rank as usize);
                        let (new_kmer, _) = left_extend_kmer(&kmer, None, sbwt, (suffix_len - k as u32) as usize);
                        assert_eq!(new_kmer.len(), suffix_len as usize);
                        encoding_type = 1;
                        new_kmer
                    } else {
                        encoding_type = 0;
                        sbwt.access_kmer(colex_rank as usize)
                    };

                    sequence.extend(kmer[(kmer.len() - (suffix_len as usize))..kmer.len()].iter());
                } else {
                    let length: usize = ((flag & 0b11111100) >> 2) as usize;

                    let mut arr1: [u8; 8] = [0; 8];
                    arr1[0..7].copy_from_slice(&bytes[0..7]);

                    let bitnucs = u64::from_ne_bytes(arr1);
                    let mut kmer: Vec<u8> = Vec::new();
                    let _ = bitnuc::from_2bit(bitnucs, length, &mut kmer);
                    bases += kmer.len();
                    sequence.extend(kmer.iter());
                    encoding_type = 3;
                    colex = None;
                }

                let end = bases;

                let contig_id = sequences.len();
                positions.push((contig_id, start, end, end - start, encoding_type, colex));
                if first {
                    bases = 0;
                    sequences.push(sequence.clone());
                    sequence.clear();
                }
            });
        },
    }

    positions
}

pub fn read_block<R: std::io::Read>(
    _file_header: &FileHeader,
    _sbwt: &SbwtIndexVariant,
    conn: &mut R,
) -> Result<Vec<u64>, E> {
    // Colex ranks
    let mut header_bytes_1: [u8; 32] = [0_u8; 32];
    conn.read_exact(&mut header_bytes_1)?;
    let header_1 = decode_block_header(&header_bytes_1)?;

    let mut bytes_1: Vec<u8> = vec![0; header_1.block_size as usize];
    let _ = conn.read_exact(&mut bytes_1);

    // Match lengths
    let mut header_bytes_2: [u8; 32] = [0; 32];
    let _ = conn.read_exact(&mut header_bytes_2);
    let header_2 = decode_block_header(&header_bytes_2)?;

    let mut bytes_2: Vec<u8> = vec![0; header_2.block_size as usize];
    let _ = conn.read_exact(&mut bytes_2);

    // Flags
    let mut header_bytes_3: [u8; 32] = [0; 32];
    let _ = conn.read_exact(&mut header_bytes_3);
    let header_3 = decode_block_header(&header_bytes_3)?;

    let mut bytes_3: Vec<u8> = vec![0; header_3.block_size as usize];
    let _ = conn.read_exact(&mut bytes_3);

    // Bitnuc coded data
    let mut header_bytes_4: [u8; 32] = [0; 32];
    let _ = conn.read_exact(&mut header_bytes_4);
    let header_4 = decode_block_header(&header_bytes_4)?;

    let mut bytes_4: Vec<u8> = vec![0; header_4.block_size as usize];
    let _ = conn.read_exact(&mut bytes_4);

    // Decompress
    let colex_ranks = decode::decompress_block(&bytes_1, &header_1, crate::encode::Codec::MinimalBinary)?;
    let match_lengths = decode::decompress_block(&bytes_2, &header_2, crate::encode::Codec::Rice)?;
    let flags = decode::decompress_block(&bytes_3, &header_3, crate::encode::Codec::Rice)?;
    let bitnuc_codings = decode::decompress_block(&bytes_4, &header_4, crate::encode::Codec::MinimalBinary)?;

    let decompressed: Vec<u64> = decode::zip_block_contents(&colex_ranks, &match_lengths, &flags, &bitnuc_codings)?;

    Ok(decompressed)
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
    fn left_extensions_with_and_without_lcs_agree() {
        // Related sequences with substitutions, so the de Bruijn graph branches.
        let base = random_sequence(3000, 1);
        let mut seqs = vec![base.clone()];
        for (seed, step) in [(2, 97), (3, 211), (4, 401)] {
            let mut variant = base.clone();
            let flips = random_sequence(variant.len() / step + 1, seed);
            for (n, i) in (step / 2..variant.len()).step_by(step).enumerate() {
                if flips[n] != variant[i] { variant[i] = flips[n]; }
            }
            seqs.push(variant);
        }
        let mut opts = kbo::BuildOpts::default();
        opts.k = 15;
        opts.add_revcomp = true;
        opts.build_select = true;
        let (variant, lcs) = kbo::build(&seqs, opts);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = &variant;
        let k = sbwt.k();

        let mut n_extended = 0;
        for seq in &seqs {
            for end in k..=seq.len() {
                let kmer = &seq[end - k..end];
                let max = end - k;
                let searched = left_extend_kmer(kmer, None, sbwt, max);
                let plain = searched.0.clone();
                // Walking the SBWT from the k-mer's row, with the LCS array or the
                // LCS values from the SBWT: checked, and by a known length.
                let row = sbwt.search(kmer).unwrap().start;
                let (with_array, without) = (walk::LcsIndex { sbwt, lcs: &lcs }, walk::WalkIndex { sbwt, lcs: None });
                assert_eq!(with_array.left_extend(row, max, true), searched);
                assert_eq!(without.left_extend(row, max, true), searched);
                assert_eq!(with_array.left_extend(row, plain.len() - k, false).0, plain);
                n_extended += (plain.len() > k) as usize;
            }

            // The encoder's extensions from the matching statistics are those of
            // left_extend_kmer with the sequence as the reference.
            let (dictionary, implicit) = encode_sequence(seq, &variant, &lcs).unwrap();
            let mut end = seq.len();
            for ((len, _), &imp) in dictionary.iter().zip(&implicit) {
                if *len >= k && end > k + 1 {
                    let extended = left_extend_kmer(&seq[end - k..end], Some(&seq[..end]), sbwt, end - k - 1);
                    assert_eq!(extended, (seq[end - len..end].to_vec(), imp));
                }
                end -= len;
            }

            // Both decoders recover the sequence from its encoding.
            let path = |_| dictionary.iter().zip(&implicit).map(|((len, colex), &imp)| {
                let len = *len as u32;
                let suffix_len = if imp {
                    graph::SuffixLen::Implicit
                } else if len > k as u32 {
                    graph::SuffixLen::Explicit { suffix_len: len }
                } else {
                    graph::SuffixLen::Plain { suffix_len: len }
                };
                (colex.start as u32, suffix_len)
            }).collect::<Vec<_>>();
            assert_eq!(&graph::decode_path(path(()), seq.len(), &variant, None).unwrap(), seq);
            assert_eq!(&graph::decode_path(path(()), seq.len(), &variant, Some(&lcs)).unwrap(), seq);
        }
        assert!(n_extended > 0, "no k-mer was left extended");
    }

    fn related_index(k: usize) -> (Vec<Vec<u8>>, SbwtIndexVariant, LcsArray) {
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
        let mut opts = kbo::BuildOpts::default();
        opts.k = k;
        opts.add_revcomp = true;
        opts.build_select = true;
        let (variant, lcs) = kbo::build(&seqs, opts);
        (seqs, variant, lcs)
    }

    #[test]
    fn lcs_array_and_lcs_from_sbwt_agree() {
        for k in [15, 31] {
            let (_, variant, lcs) = related_index(k);
            let SbwtIndexVariant::SubsetMatrix(sbwt) = &variant;
            assert_eq!(lcs.len(), sbwt.n_sets());
            for rank in 0..lcs.len() {
                assert_eq!(walk::lcs_from_sbwt(sbwt, rank), lcs.access(rank), "k {} rank {}", k, rank);
            }
        }
    }

    #[test]
    fn lcs_context_suffix_codes_roundtrip() {
        let (seqs, variant, lcs) = related_index(15);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = &variant;
        let mut builder = graph::GraphBuilder::new();
        for seq in &seqs {
            let (dictionary, implicit) = encode_sequence(seq, &variant, &lcs).unwrap();
            builder.add_path_with_flags(seq.len(), &dictionary, &implicit);
        }
        let (chains, _) = builder.to_chains();
        let from_array = |rank: u32| lcs.access(rank as usize) as u32;
        let from_sbwt = |rank: u32| walk::lcs_from_sbwt(sbwt, rank as usize) as u32;
        for cap in [1, 8, 255] {
            let stored = chains.to_stored_with(Some(&from_array), cap, Some(15), None);
            let bytes = postcard::to_allocvec(&stored).unwrap();
            // Encoded with the array, decoded with values computed from the SBWT (and vice versa).
            let decoded = graph::collapsed::ChainCsr::from_stored_with(postcard::from_bytes(&bytes).unwrap(), Some(&from_sbwt), Some(15), None).unwrap();
            assert_eq!(decoded, chains, "cap {}", cap);
            let stored = chains.to_stored_with(Some(&from_sbwt), cap, Some(15), None);
            assert_eq!(graph::collapsed::ChainCsr::from_stored_with(stored, Some(&from_array), Some(15), None).unwrap(), chains);
            // The LCS values are required to decode.
            let stored: graph::collapsed::StoredChains = postcard::from_bytes(&bytes).unwrap();
            assert!(graph::collapsed::ChainCsr::from_stored_with(stored, None, Some(15), None).is_err());
        }
        // Without contexts nothing changes.
        assert_eq!(graph::collapsed::ChainCsr::from_stored(chains.to_stored()).unwrap(), chains);
    }

    #[test]
    fn nodes_derived_by_walks_roundtrip() {
        use graph::GraphWalk;
        let (seqs, variant, lcs) = related_index(15);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = &variant;
        let mut builder = graph::GraphBuilder::new();
        for seq in &seqs {
            let (dictionary, implicit) = encode_sequence(seq, &variant, &lcs).unwrap();
            builder.add_path_with_flags(seq.len(), &dictionary, &implicit);
        }
        let (chains, _) = builder.to_chains();
        let from_array = |rank: u32| lcs.access(rank as usize) as u32;
        let from_sbwt = |rank: u32| walk::lcs_from_sbwt(sbwt, rank as usize) as u32;
        let walk = graph::collapsed::WalkIndex { sbwt, lcs: Some(&lcs) };
        let plain = chains.to_stored_with(Some(&from_array), 255, Some(15), None);
        let stored = chains.to_stored_with(Some(&from_array), 255, Some(15), Some(walk));
        assert!(stored.derived && !plain.derived);
        // Some nodes are derived: fewer ranks are stored.
        assert!(stored.node_ranks.len() < plain.node_ranks.len(), "{} vs {}", stored.node_ranks.len(), plain.node_ranks.len());
        // With walks the targets are coded where each chain ends, some of
        // them junctions reached before.
        assert!(stored.targets.is_empty() && !plain.targets.is_empty());
        let mut reached = vec![false; chains.junction_weights.len()];
        let mut reached_before = 0;
        for junction in 0..chains.junction_weights.len() {
            let (start, end) = chains.row(junction);
            for chain in start..end {
                let target = chains.targets[chain] as usize;
                reached_before += reached[target] as usize;
                reached[target] = true;
            }
        }
        assert!(reached_before > 0, "no chain ends at a junction reached before");
        let bytes = postcard::to_allocvec(&stored).unwrap();
        // Decoded with the LCS array, and with LCS values computed from the SBWT (no array).
        let decoded = graph::collapsed::ChainCsr::from_stored_with(postcard::from_bytes(&bytes).unwrap(), Some(&from_array), Some(15), Some(walk)).unwrap();
        assert_eq!(decoded, chains);
        let no_array = graph::collapsed::WalkIndex { sbwt, lcs: None };
        let decoded = graph::collapsed::ChainCsr::from_stored_with(postcard::from_bytes(&bytes).unwrap(), Some(&from_sbwt), Some(15), Some(no_array)).unwrap();
        assert_eq!(decoded, chains);
        // The derived nodes need the index.
        assert!(graph::collapsed::ChainCsr::from_stored_with(postcard::from_bytes(&bytes).unwrap(), Some(&from_array), Some(15), None).is_err());
        // A corrupt walk stream is an error, not a panic.
        let mut bad: graph::collapsed::StoredChains = postcard::from_bytes(&bytes).unwrap();
        bad.node_walks.truncate(bad.node_walks.len() / 2);
        let _ = graph::collapsed::ChainCsr::from_stored_with(bad, Some(&from_array), Some(15), Some(walk));
    }

    #[test]
    fn archives_written_with_the_lcs_array_need_an_lcs_source() {
        use graph::{LcsSource, StoredGraph, GraphWalk};
        let (seqs, variant, lcs) = related_index(15);
        let mut builder = graph::GraphBuilder::new().collapse_chains(true);
        for seq in &seqs {
            builder.add_contig(seq, &variant, &lcs).unwrap();
        }
        let (stored, choices) = builder.finish();
        let walks: Vec<_> = (0..seqs.len() as u32).map(|c| stored.walk_path(&choices, c).unwrap()).collect();
        let mut with_lcs = Vec::new();
        graph::encoder::GraphEncoder::new(stored, choices).with_index(&variant, &lcs).encode_to_write(&mut with_lcs).unwrap();

        for source in [LcsSource::Array(&lcs), LcsSource::Sbwt(&variant)] {
            let mut input = with_lcs.as_slice();
            let (graph, choices, ..) = graph::decoder::GraphDecoder::new(&mut input).unwrap().decode_from_read_with(Some(&variant), Some(source), Some(15)).unwrap();
            for (color, walk) in walks.iter().enumerate() {
                assert_eq!(&graph.walk_path(&choices, color as u32).unwrap(), walk);
            }
            assert!(matches!(graph, StoredGraph::Chains(_)));
        }
        let mut input = with_lcs.as_slice();
        assert!(graph::decoder::GraphDecoder::new(&mut input).unwrap().decode_from_read().is_err());
    }
}
