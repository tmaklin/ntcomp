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
use sbwt::SbwtIndexVariant;

pub mod encode;
pub mod decode;


#[derive(Encode, Decode)]
pub struct HeaderPlaceholder {
    pub ph1: u64,
    pub ph2: u64,
    pub ph3: u64,
    pub ph4: u64,
}

#[derive(Encode, Decode)]
pub struct BlockHeader {
    pub block_size: u32,
    pub num_records: u32,
    pub num_u64: u32,
    pub encoded_size: u32,
    pub rice_param: u8,
    pub bitpacker_exponent: u8,
    pub placeholder1: u64,
    pub placeholder2: u32,
    pub placeholder3: u16,
    // TODO should store the used codec so others can be used
    // TODO store total number of bases?
}

pub fn encode_file_header(
    ph1: u64,
    ph2: u64,
    ph3: u64,
    ph4: u64
) -> Vec<u8> {
    let mut bytes: Vec<u8> = Vec::new();
    let header_placeholder = HeaderPlaceholder{ ph1, ph2, ph3, ph4 };
    let nbytes = encode_into_std_write(
        &header_placeholder,
        &mut bytes,
        bincode::config::standard().with_fixed_int_encoding(),
    );
    assert_eq!(nbytes.unwrap(), 32);
    bytes
}

pub fn decode_file_header(
    header_bytes: &[u8],
) -> HeaderPlaceholder {
    decode_from_slice(header_bytes, bincode::config::standard().with_fixed_int_encoding()).unwrap().0
}

pub fn encode_block_header(
    header: &BlockHeader,
) -> Vec<u8> {
    let mut bytes: Vec<u8> = Vec::new();
    let nbytes = encode_into_std_write(
        header,
        &mut bytes,
        bincode::config::standard().with_fixed_int_encoding(),
    );
    assert_eq!(nbytes.unwrap(), 32);
    bytes
}

pub fn decode_block_header(
    header_bytes: &[u8],
) -> BlockHeader {
    decode_from_slice(header_bytes, bincode::config::standard().with_fixed_int_encoding()).unwrap().0
}

pub fn left_extend_kmer(
    kmer_start: &[u8],
    ref_nts: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    max_extension_len: usize,
) -> Vec<u8> {
    assert!(!kmer_start.is_empty());
    assert!(!ref_nts.is_empty());

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
            let seq_matches = new_kmers[0].0[0] == ref_nts[ref_nts.len() - kmer.len() - 1];
            if seq_matches && new_kmers.len() == 1 && new_kmers[0].1.end - new_kmers[0].1.start == 1 {
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

pub fn left_extend_kmer2(
    kmer_start: &[u8],
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    max_extension_len: usize,
) -> Vec<u8> {
    assert!(!kmer_start.is_empty());

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

pub fn encode_sequence(
    nucleotides: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> Vec<u64> {
    let n = nucleotides.len();
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
                    let new_kmer = left_extend_kmer(&kmer, ref_nts, sbwt, i - k - 1);
                    let new_seq = nucleotides[(i - new_kmer.len())..i].to_vec();
                    assert_eq!(new_kmer.len(), new_seq.len());

                    let mut match_len = new_kmer.len();

                    let old_i = i;
                    loop {
                        if res[i - 1].0 < match_len {
                            let old_ms = res[i - 1].0;
                            i -= old_ms;
                            match_len -= old_ms;
                        } else {
                            res[old_i - 1].0 = new_kmer.len() - (match_len - 1);
                            kept.push(old_i - 1);
                            break;
                        }
                    }
                } else {
                    kept.push(i - 1);
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
    let dictionary_max = dictionary.iter().map(|x| x.0).max().unwrap();

    assert!(dictionary_max < 16777216); // Match lengths are coded as 24 bit unsigned integers
    assert_eq!(dictionary_bases, nucleotides.len());

    encode::encode_dictionary(&dictionary)
}

pub fn write_block_to<W: std::io::Write>(
    u64_encoding: &[u64],
    num_records: usize,
    sink: &mut W,
) -> Result<(), Box<dyn std::error::Error>> {
    let (data_1, data_2) = encode::split_encoded_dictionary(u64_encoding);

    let block_1 = encode::compress_block(&data_1, num_records, encode::Codec::Rice);
    let block_2 = encode::compress_block(&data_2, num_records, encode::Codec::Rice);

    sink.write_all(&block_1)?;
    sink.write_all(&block_2)?;
    Ok(())
}

pub fn decode_sequence(
    dictionary: &[(u32, u32, bool)],
    sbwt: &SbwtIndexVariant,
) -> Vec<u8> {
    let mut sequence: Vec<u8> = Vec::new();
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            dictionary.iter().for_each(|record| {
                let (suffix_len, colex_rank, _) = *record;
                let kmer = if suffix_len > k as u32 {
                    let kmer = sbwt.access_kmer(colex_rank as usize);
                    let new_kmer = left_extend_kmer2(&kmer, sbwt, (suffix_len - k as u32) as usize);
                    assert_eq!(new_kmer.len(), suffix_len as usize);
                    new_kmer
                } else {
                    sbwt.access_kmer(colex_rank as usize)
                };
                sequence.extend(kmer[(kmer.len() - (suffix_len as usize))..kmer.len()].iter().rev());
            });
        },
    }
    sequence.into_iter().rev().collect()
}
