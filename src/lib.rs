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

use kbo::variant_calling::Variant;

use sbwt::LcsArray;
use sbwt::StreamingIndex;
use sbwt::SbwtIndexVariant;

pub mod encode;
pub mod decode;
pub mod delta_encode;

type E = Box<dyn std::error::Error>;

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
    pub rice_param: u64,
    pub bitpacker_exponent: u8,
    pub placeholder1: u8,
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
) -> Result<Vec<u8>, E> {
    let mut bytes: Vec<u8> = Vec::new();
    let header_placeholder = HeaderPlaceholder{ ph1, ph2, ph3, ph4 };
    let nbytes = encode_into_std_write(
        &header_placeholder,
        &mut bytes,
        bincode::config::standard().with_fixed_int_encoding(),
    )?;
    assert_eq!(nbytes, 32);
    Ok(bytes)
}

pub fn decode_file_header(
    header_bytes: &[u8],
) -> Result<HeaderPlaceholder, E> {
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
) -> Result<Vec<(usize, Range<usize>)>, E> {
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
    let dictionary_max = dictionary.iter().map(|x| x.0).max().ok_or(encode::EncodeError{})?;

    assert!(dictionary_max < 16777216); // Match lengths are coded as 24 bit unsigned integers
    assert_eq!(dictionary_bases, nucleotides.len());

    Ok(dictionary)
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
                        let new_kmer = left_extend_kmer2(&kmer, sbwt, (suffix_len - k as u32) as usize);
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

pub fn decode_block<R: std::io::Read>(
    _file_header: &HeaderPlaceholder,
    sbwt: &SbwtIndexVariant,
    conn: &mut R,
) -> Result<Vec<Vec<u8>>, E> {
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

    // Variants
    // Positions
    let mut header_bytes_5: [u8; 32] = [0_u8; 32];
    conn.read_exact(&mut header_bytes_5)?;
    let header_5 = decode_block_header(&header_bytes_5)?;

    let mut bytes_5: Vec<u8> = vec![0; header_5.block_size as usize];
    let _ = conn.read_exact(&mut bytes_5);

    // Query chars
    let mut header_bytes_6: [u8; 32] = [0; 32];
    let _ = conn.read_exact(&mut header_bytes_6);
    let header_6 = decode_block_header(&header_bytes_6)?;

    let mut bytes_6: Vec<u8> = vec![0; header_6.block_size as usize];
    let _ = conn.read_exact(&mut bytes_6);

    // Ref variant lengths
    let mut header_bytes_7: [u8; 32] = [0; 32];
    let _ = conn.read_exact(&mut header_bytes_7);
    let header_7 = decode_block_header(&header_bytes_7)?;

    let mut bytes_7: Vec<u8> = vec![0; header_7.block_size as usize];
    let _ = conn.read_exact(&mut bytes_7);

    // Query variant lengths
    let mut header_bytes_8: [u8; 32] = [0; 32];
    let _ = conn.read_exact(&mut header_bytes_8);
    let header_8 = decode_block_header(&header_bytes_8)?;

    let mut bytes_8: Vec<u8> = vec![0; header_8.block_size as usize];
    let _ = conn.read_exact(&mut bytes_8);

    // Num variants per seq record
    let mut header_bytes_9: [u8; 32] = [0; 32];
    let _ = conn.read_exact(&mut header_bytes_9);
    let header_9 = decode_block_header(&header_bytes_9)?;

    let mut bytes_9: Vec<u8> = vec![0; header_9.block_size as usize];
    let _ = conn.read_exact(&mut bytes_9);


    // Decompress
    let decompressed_1 = decode::decompress_block(&bytes_1, &header_1, crate::encode::Codec::MinimalBinary)?;
    let decompressed_2 = decode::decompress_block(&bytes_2, &header_2, crate::encode::Codec::Rice)?;
    let decompressed_3 = decode::decompress_block(&bytes_3, &header_3, crate::encode::Codec::Rice)?;
    let decompressed_4 = decode::decompress_block(&bytes_4, &header_4, crate::encode::Codec::MinimalBinary)?;

    let decompressed_5 = decode::decompress_block(&bytes_5, &header_5, crate::encode::Codec::Rice)?;
    let decompressed_6 = decode::decompress_block(&bytes_6, &header_6, crate::encode::Codec::MinimalBinary)?;
    let decompressed_7= decode::decompress_block(&bytes_7, &header_7, crate::encode::Codec::Rice)?;
    let decompressed_8 = decode::decompress_block(&bytes_8, &header_8, crate::encode::Codec::Rice)?;
    let decompressed_9 = decode::decompress_block(&bytes_9, &header_9, crate::encode::Codec::Rice)?;

    let decompressed: Vec<u64> = decode::zip_block_contents(&decompressed_1, &decompressed_2, &decompressed_3, &decompressed_4)?;
    let variants: Vec<Variant> = delta_encode::decode(&decompressed_5, &decompressed_6, &decompressed_7, &decompressed_8);

    let mut decoded = decode_sequence(&decompressed, sbwt);

    assert_eq!(decoded.len(), decompressed_9.len());

    let mut i: usize = 0;
    decoded = decompressed_9.iter().enumerate().map(|(idx, num_variants)| {
        let seq_variants = variants[i..(i + *num_variants as usize)].to_vec();
        let seq = kbo::variant_calling::revert_edits(&decoded[idx], &seq_variants);
        i += *num_variants as usize;
        seq
    }).collect();

    Ok(decoded)
}
