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
    // The different lenght variants take a lot of space for some reason
    let variants: Vec<Variant> = variants.iter().filter_map(|x| {
        if x.ref_chars.is_empty() || x.query_chars.is_empty() || x.ref_chars.len() == x.query_chars.len() {
            Some(x.clone())
        } else {
            None
        }
    }).collect();

    let positions: Vec<u64> = variants.iter().map(|x| x.query_pos as u64).collect();
    let q_nts: Vec<u8> = variants.iter().flat_map(|x| x.query_chars.clone()).collect();
    let r_lens: Vec<u64> = variants.iter().map(|x| x.ref_chars.len() as u64).collect();

    let enc1 = crate::encode::rice_encode(&positions).unwrap().0;
    let q_bitnucs: Vec<u64> = q_nts.chunks(31).map(|x| bitnuc::as_2bit(x).unwrap()).collect();
    let enc2 = crate::encode::minimal_binary_encode(&q_bitnucs).unwrap().0;
    let enc3 = crate::encode::rice_encode(&r_lens).unwrap().0;

    let mut block1: Vec<u8> = crate::encode::deflate_bytes(&enc1.iter().flat_map(|x| {
        let mut arr: [u8; 8] = [0; 8];
        arr.copy_from_slice(&x.to_ne_bytes()[0..8]);
        arr
    }).collect::<Vec<u8>>()).unwrap();

    let block2: Vec<u8> = crate::encode::deflate_bytes(&enc2.iter().flat_map(|x| {
        let mut arr: [u8; 8] = [0; 8];
        arr.copy_from_slice(&x.to_ne_bytes()[0..8]);
        arr
    }).collect::<Vec<u8>>()).unwrap();

    let block3: Vec<u8> = crate::encode::deflate_bytes(&enc3.iter().flat_map(|x| {
        let mut arr: [u8; 8] = [0; 8];
        arr.copy_from_slice(&x.to_ne_bytes()[0..8]);
        arr
    }).collect::<Vec<u8>>()).unwrap();

    block1.extend(block2.iter());
    block1.extend(block3.iter());
    block1
}
