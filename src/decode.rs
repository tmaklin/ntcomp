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
use std::io::Write;

use dsi_bitstream::traits::BE;
use dsi_bitstream::prelude::MemWordReader;
use dsi_bitstream::prelude::BufBitReader;
use dsi_bitstream::codes::RiceRead;
use flate2::write::GzDecoder;

use crate::encode;

pub fn decode_block(
    block: &[u8],
    header: &encode::BlockHeader,
) -> Vec<(u8, u32, bool)> {
    // info!("Inflating block...");
    let mut inflated: Vec<u8> = Vec::new();
    let mut decoder = GzDecoder::new(&mut inflated);
    decoder.write_all(block).unwrap();
    inflated = decoder.finish().unwrap().to_vec();
    let bytes = inflated;

    // info!("Converting u32 to u64...");
    let rice_encoded: Vec<u64> = bytes.chunks(8).map(|x| {
        // Convert bytes to Golomb-Rice encoded u64
        let arr: [u8; 8] = [x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]];
        // Decode
        u64::from_ne_bytes(arr)
    }).collect();
    assert_eq!(rice_encoded.len(), header.encoded_size as usize);

    // info!("Decoding RiceCodec...");
    let mut reader = BufBitReader::<BE, _>::new(MemWordReader::new(rice_encoded));
    let encoded: Vec<u64> = (0..header.num_u64).map(|_| {
        reader.read_rice(header.rice_param as usize).unwrap()
    }).collect();

    // info!("Decoding u64 encoded (MS, colex interval) pairs...");
    let decoded: Vec<(u8, u32, bool)> = encoded.iter().map(|x| {
        let arr: [u8; 8] = x.to_ne_bytes();
        let mut arr1: [u8; 4] = [0; 4];
        let mut arr2: [u8; 1] = [0; 1];
        let mut arr3: [u8; 1] = [0; 1];
        arr1.copy_from_slice(&arr[0..4]);
        arr2.copy_from_slice(&arr[4..5]);
        arr3.copy_from_slice(&arr[5..6]);
        let colex_rank = u32::from_ne_bytes(arr1);
        let ms = u8::from_ne_bytes(arr2);
        (ms, colex_rank, u8::from_ne_bytes(arr3) == 1)
    }).collect();

    decoded
}
