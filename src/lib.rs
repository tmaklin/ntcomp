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

use sbwt::LcsArray;
use sbwt::StreamingIndex;
use sbwt::SbwtIndex;
use sbwt::SbwtIndexVariant;
use sbwt::SubsetMatrix;

pub mod encode;
pub mod decode;

pub fn encode_sequence(
    nucleotides: &[u8],
    index: &StreamingIndex<SbwtIndex<SubsetMatrix>, LcsArray>,
) -> Vec<u64> {
    let dictionary: Vec<(usize, Range<usize>)> = index.matching_statistics(nucleotides);
    encode::encode_dictionary(&dictionary)
}

pub fn decode_sequence(
    dictionary: &[(u8, u32, bool)],
    sbwt: &SbwtIndexVariant,
//    sbwt: &SbwtIndex<SubsetMatrix>,
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
