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
//! The plain NLZ file format: the dictionaries ([super::parse]) of the
//! sequences one after another, in blocks of records split into four
//! streams (colex ranks, match lengths, flags, 2-bit coded short entries),
//! each Rice or minimal binary coded and deflated ([encode], [decode]).
//!
//! A file starts with [HEADER_BYTES] zero bytes (a header that was never
//! filled in), then the blocks. Decoding extends the _k_-mers with
//! [left_extend_kmer] (searches in the SBWT).
//!
use bincode::{Encode, Decode};
use bincode::encode_into_std_write;
use bincode::decode_from_slice;

use sbwt::sbwt_index_variant::SbwtIndexVariant;

use super::left_extend_kmer;

pub mod decode;
pub mod encode;

type E = Box<dyn std::error::Error>;

/// Size of the (unused, all zero) file header.
pub const HEADER_BYTES: usize = 50;

/// Writes the file header.
pub fn write_header<W: std::io::Write>(
    sink: &mut W,
) -> Result<(), E> {
    sink.write_all(&[0_u8; HEADER_BYTES])?;
    Ok(())
}

/// Reads past the file header.
pub fn read_header<R: std::io::Read>(
    conn: &mut R,
) -> Result<(), E> {
    let mut header = [0_u8; HEADER_BYTES];
    conn.read_exact(&mut header)?;
    Ok(())
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

/// Writes the encoded dictionaries of `num_records` sequences as a block.
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

/// Decodes the sequences of a block ([read_block]).
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

/// Reads the next block written by [write_block_to].
pub fn read_block<R: std::io::Read>(
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
    let colex_ranks = decode::decompress_block(&bytes_1, &header_1, encode::Codec::MinimalBinary)?;
    let match_lengths = decode::decompress_block(&bytes_2, &header_2, encode::Codec::Rice)?;
    let flags = decode::decompress_block(&bytes_3, &header_3, encode::Codec::Rice)?;
    let bitnuc_codings = decode::decompress_block(&bytes_4, &header_4, encode::Codec::MinimalBinary)?;

    let decompressed: Vec<u64> = decode::zip_block_contents(&colex_ranks, &match_lengths, &flags, &bitnuc_codings)?;

    Ok(decompressed)
}
