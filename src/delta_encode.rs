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
use crate::encode::Codec;
use crate::encode::compress_block;

use kbo::variant_calling::Variant;

use sbwt::LcsArray;
use sbwt::SbwtIndexVariant;

pub fn compute_deltas(
    seqrec: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> (Vec<u8>, Vec<Variant>) {
    let mut call_opts =  kbo::CallOpts::default();
    call_opts.sbwt_build_opts.k = 255;

    let variants = kbo::call(sbwt, lcs, seqrec, call_opts.clone());

    let edited = kbo::variant_calling::edit(seqrec, &variants);

    (edited, variants)
}

pub fn encode(
    variants_v: &[Vec<Variant>],
) -> Vec<u8> {
    let variants: Vec<Variant> = variants_v.iter().cloned().flatten().collect();
    let positions: Vec<u64> = variants.iter().map(|x| x.query_pos as u64).collect();
    let q_nts: Vec<u8> = variants.iter().flat_map(|x| x.query_chars.clone()).collect();
    let r_lens: Vec<u64> = variants.iter().map(|x| x.ref_chars.len() as u64).collect();
    let q_lens: Vec<u64> = variants.iter().map(|x| x.query_chars.len() as u64).collect();
    let num_records: Vec<u64> = variants_v.iter().map(|x| x.len() as u64).collect();

    let q_bitnucs: Vec<u64> = q_nts.chunks(31).map(|x| bitnuc::as_2bit(x).unwrap()).collect();

    let mut blocks = compress_block(&positions, variants.len(), Codec::Rice).unwrap();
    blocks.extend(compress_block(&q_bitnucs, variants.len(), Codec::MinimalBinary).unwrap().iter());
    blocks.extend(compress_block(&r_lens, variants.len(), Codec::Rice).unwrap().iter());
    blocks.extend(compress_block(&q_lens, variants.len(), Codec::Rice).unwrap().iter());
    blocks.extend(compress_block(&num_records, num_records.len(), Codec::Rice).unwrap().iter());

    blocks
}

pub fn decode(
    positions: &[u64],
    q_bitnucs: &[u64],
    r_lens: &[u64],
    q_lens: &[u64],
) -> Vec<Variant> {
    let bitnuc_len = q_lens.iter().sum::<u64>();
    let bitnuc_vals: Vec<u8> = q_bitnucs.iter().enumerate().flat_map(|(idx, bitnucs)| {
        let len = if idx + 1 == q_bitnucs.len() {
            bitnuc_len % 31
        } else {
            31
        };
        let mut kmer: Vec<u8> = Vec::new();
        bitnuc::from_2bit(*bitnucs, len as usize, &mut kmer).unwrap();
        kmer
    }).collect();

    let mut j = 0;
    let res: Vec<Variant> = q_lens.iter().enumerate().map(|(idx, length)| {
        let q_nts: Vec<u8> = bitnuc_vals[j..(j + *length as usize)].to_vec();
        j += *length as usize;
        Variant{ query_chars: q_nts, ref_chars: vec![0_u8; r_lens[idx] as usize], query_pos: positions[idx] as usize }
    }).collect::<Vec<Variant>>();

    res
}
