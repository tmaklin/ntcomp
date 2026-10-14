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

pub mod encode;
pub mod decode;
pub mod graph;

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

/// Left extends `kmer_start` while the extension is unique in the SBWT and
/// matches `ref_nts`, by at most `max_extension_len` characters.
///
/// Returns the extended k-mer and whether the extension stopped because the
/// SBWT has no unique left extension, i.e. where [left_extend_kmer2] stops as
/// well. False if it stopped because of a mismatch with `ref_nts` or because
/// `max_extension_len` was reached.
///
/// If `lcs` is given, the extension is computed as right extensions of the
/// reverse complement, see [left_extend_rc]: one search plus constant time per
/// character instead of four searches per character. This requires that the
/// SBWT contains the reverse complements of its _k_-mers (`ntcomp build` adds
/// them), and gives the same result as without `lcs`.
pub fn left_extend_kmer(
    kmer_start: &[u8],
    ref_nts: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: Option<&LcsArray>,
    max_extension_len: usize,
) -> (Vec<u8>, bool) {
    assert!(!kmer_start.is_empty());
    assert!(!ref_nts.is_empty());

    if let Some(extended) = lcs.and_then(|lcs| left_extend_rc(kmer_start, Some(ref_nts), sbwt, lcs, max_extension_len)) {
        return extended;
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
            let seq_matches = new_kmers[0].0[0] == ref_nts[ref_nts.len() - kmer.len() - 1];
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

/// Left extends `kmer_start` while the extension is unique in the SBWT, by at
/// most `max_extension_len` characters.
///
/// If `lcs` is given, uses right extensions of the reverse complement like
/// [left_extend_kmer].
pub fn left_extend_kmer2(
    kmer_start: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: Option<&LcsArray>,
    max_extension_len: usize,
) -> Vec<u8> {
    assert!(!kmer_start.is_empty());

    if let Some((extended, _)) = lcs.and_then(|lcs| left_extend_rc(kmer_start, None, sbwt, lcs, max_extension_len)) {
        return extended;
    }

    let mut left_extension_len = 0;
    let mut kmer = kmer_start.to_vec();
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
            if new_kmers.len() == 1 && new_kmers[0].1.end - new_kmers[0].1.start == 1 {
                kmer = [&[new_kmers[0].0[0]], kmer.as_slice()].concat();
            } else {
                break;
            }
        } else {
            break;
        }
        left_extension_len += 1;
    }
    kmer
}

/// Left extension of `kmer_start` as right extensions of its reverse complement.
///
/// The predecessors of a _k_-mer are the reverse complements of the successors
/// of its reverse complement if the SBWT contains both strands. Finds the
/// reverse complement with one search, then takes one
/// [extend_right](sbwt::ExtendRight::extend_right) per character and base.
/// Only the first _k_-mer in colex order among those sharing the same _k_-1
/// suffix stores the outgoing edges in the SBWT, so the extensions start from
/// that group, found with the `lcs` array.
///
/// Stops under the same conditions as [left_extend_kmer] (with `ref_nts`) or
/// [left_extend_kmer2] (without) and returns the same values. Returns None if
/// the reverse complement of `kmer_start` is not in the SBWT.
fn left_extend_rc(
    kmer_start: &[u8],
    ref_nts: Option<&[u8]>,
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: &LcsArray,
    max_extension_len: usize,
) -> Option<(Vec<u8>, bool)> {
    use sbwt::ExtendRight;

    let k = sbwt.k();
    let mut extension: Vec<u8> = Vec::new();
    let mut stopped_in_sbwt = false;
    if max_extension_len > 0 {
        let mut colex = sbwt.search(&graph::orientation::reverse_complement(kmer_start))?.start;
        while extension.len() < max_extension_len {
            let mut group_start = colex;
            while group_start > 0 && lcs.access(group_start) >= k - 1 {
                group_start -= 1;
            }
            let mut successor: Option<(u8, usize)> = None;
            let mut n_successors = 0;
            for &c in sbwt.alphabet() {
                let res = sbwt.extend_right(group_start..(colex + 1), c);
                if !res.is_empty() {
                    n_successors += 1;
                    successor = Some((c, res.start));
                }
            }
            if n_successors != 1 {
                stopped_in_sbwt = true;
                break;
            }
            let (c, next) = successor.unwrap();
            let predecessor = match c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            };
            if ref_nts.is_some_and(|ref_nts| predecessor != ref_nts[ref_nts.len() - kmer_start.len() - extension.len() - 1]) {
                break;
            }
            extension.push(predecessor);
            colex = next;
        }
    }
    extension.reverse();
    extension.extend_from_slice(kmer_start);
    Some((extension, stopped_in_sbwt))
}

pub fn encode_sequence(
    nucleotides: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> Result<Vec<(usize, Range<usize>)>, E> {
    Ok(encode_sequence_flagged(nucleotides, sbwt, lcs)?.0)
}

/// Same as [encode_sequence], but also returns for every dictionary entry
/// whether its suffix length is *implicit*: the entry is a left extended k-mer
/// whose extension stopped where [left_extend_kmer2] stops as well, so the
/// decoder can recover the length by extending the k-mer as far as possible.
pub fn encode_sequence_flagged(
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

            let mut i = n;
            let mut kept: Vec<usize> = Vec::new();
            while i > 0 {
                let colex_int = res[i - 1].1.clone();

                if res[i - 1].0 == k && i > k + 1 {
                    let kmer = sbwt.access_kmer(colex_int.start);
                    let seq = nucleotides[(i - kmer.len())..i].to_vec();
                    assert_eq!(kmer, seq);

                    let ref_nts = &nucleotides[0..i];
                    let (new_kmer, stopped_in_sbwt) = left_extend_kmer(&kmer, ref_nts, sbwt, Some(lcs), i - k - 1);
                    assert_eq!(new_kmer.as_slice(), &nucleotides[(i - new_kmer.len())..i]);

                    // Store the whole left extension; the next entry ends where it starts.
                    res[i - 1].0 = new_kmer.len();
                    kept.push(i - 1);
                    implicit.push(stopped_in_sbwt);
                    i -= new_kmer.len() - 1;
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
    let dictionary_max = dictionary.iter().map(|x| x.0).max().ok_or(encode::EncodeError{})?;

    assert!(dictionary_max < 16777216); // Match lengths are coded as 24 bit unsigned integers
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
                        let new_kmer = left_extend_kmer2(&kmer, sbwt, None, (suffix_len - k as u32) as usize);
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
                        let new_kmer = left_extend_kmer2(&kmer, sbwt, None, (suffix_len - k as u32) as usize);
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

