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

pub fn encode_sequence(
    nucleotides: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> Vec<u64> {
    let dictionary: Vec<(usize, Range<usize>)> = match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let index = StreamingIndex::new(sbwt, lcs);
            index.matching_statistics(nucleotides)
        },
    };
    encode::encode_dictionary(&dictionary)
}

pub fn decode_sequence(
    dictionary: &[(u8, u32, bool)],
    sbwt: &SbwtIndexVariant,
    mut pointer: usize,
) -> (Vec<u8>, usize) {
    assert!(dictionary[pointer].2);

    let mut sequence: Vec<u8> = Vec::new();
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            loop {
                let (suffix_len, colex_rank, _) = dictionary[pointer];
                let kmer = sbwt.access_kmer(colex_rank as usize);
                sequence.extend(kmer[(k - (suffix_len as usize))..k].iter().rev());
                pointer += 1;
                if pointer == dictionary.len() || dictionary[pointer].2 {
                    break;
                }
            }
        },
    }
    (sequence, pointer)
}
