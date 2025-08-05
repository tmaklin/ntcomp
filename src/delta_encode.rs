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
    variants: &[Variant],
) -> Vec<u8> {
    let positions: Vec<u64> = variants.iter().map(|x| x.query_pos as u64).collect();
    let q_nts: Vec<u8> = variants.iter().flat_map(|x| x.query_chars.clone()).collect();
    let r_lens: Vec<u64> = variants.iter().map(|x| x.ref_chars.len() as u64).collect();
    let q_lens: Vec<u64> = variants.iter().map(|x| x.query_chars.len() as u64).collect();

    let q_bitnucs: Vec<u64> = q_nts.chunks(31).map(|x| bitnuc::as_2bit(x).unwrap()).collect();

    let mut blocks = compress_block(&positions, variants.len(), Codec::Rice).unwrap();
    blocks.extend(compress_block(&q_bitnucs, variants.len(), Codec::MinimalBinary).unwrap().iter());
    blocks.extend(compress_block(&r_lens, variants.len(), Codec::Rice).unwrap().iter());
    blocks.extend(compress_block(&q_lens, variants.len(), Codec::Rice).unwrap().iter());

    blocks
}


    block1.extend(block2.iter());
    block1.extend(block3.iter());
    block1
}
