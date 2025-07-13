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

use flate2::write::GzEncoder;
use flate2::Compression;

use dsi_bitstream::impls::MemWordWriterVec;
use dsi_bitstream::impls::BufBitWriter;
use dsi_bitstream::traits::BE;
use dsi_bitstream::codes::RiceWrite;
use dsi_bitstream::codes::MinimalBinaryWrite;

use crate::BlockHeader;
use crate::encode_block_header;

fn deflate_bytes(
    bytes: &[u8],
) -> Vec<u8> {
    let mut deflated: Vec<u8> = Vec::with_capacity(bytes.len());
    let mut encoder = GzEncoder::new(&mut deflated, Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap();
    deflated
}

fn rice_encode(
    ints: &[u64],
) -> (Vec<u64>, usize) {
    let inv_mean: f64 = ((ints.len() as f64).ln() - (ints.iter().sum::<u64>() as f64).ln()).exp();
    let param: usize = dsi_bitstream::codes::rice::log2_b(inv_mean);

    let word_write = MemWordWriterVec::new(Vec::<u64>::new());
    let mut writer = BufBitWriter::<BE, _>::new(word_write);

    ints.iter().for_each(|n| { writer.write_rice(*n, param).unwrap(); } );
    let _ = writer.flush();

    (writer.into_inner().unwrap().into_inner(), param)
}

fn minimal_binary_encode(
    ints: &[u64],
) -> (Vec<u64>, usize) {
    let param = ints.iter().max().unwrap();

    let word_write = MemWordWriterVec::new(Vec::<u64>::new());
    let mut writer = BufBitWriter::<BE, _>::new(word_write);

    ints.iter().for_each(|n| { writer.write_minimal_binary(*n, *param).unwrap(); } );
    let _ = writer.flush();

    (writer.into_inner().unwrap().into_inner(), *param as usize)
}

pub fn compress_block(
    u64_encoding: &[u64],
    num_records: usize,
    codec: bool,
) -> Vec<u8> {

    let (rice_encoded, param) = if codec {
        rice_encode(u64_encoding)
    } else {
        minimal_binary_encode(u64_encoding)
    };

    let bytes = rice_encoded.iter().flat_map(|x| {
        x.to_ne_bytes()
    }).collect::<Vec<u8>>();

    let deflated: Vec<u8> = deflate_bytes(&bytes);

    let block_header = BlockHeader{ block_size: deflated.len() as u32,
                                    num_records: num_records as u32,
                                    num_u64: u64_encoding.len() as u32,
                                    encoded_size: rice_encoded.len() as u32,
                                    rice_param: param as u8,
                                    bitpacker_exponent: 8_u8,
                                    placeholder1: 0, placeholder2: 0, placeholder3: 0,
    };

    let mut block: Vec<u8> = encode_block_header(&block_header);
    block.extend(deflated.iter());
    block.shrink_to_fit();

    block
}

pub fn encode_dictionary(
    dictionary: &[(usize, std::ops::Range<usize>)],
) -> Vec<u64> {
    // TODO Should preprocess the dictionary so that the overlapping elements are removed
    let mut first: bool = true;
    let mut bases = 0;
    let mut u64_encoding: Vec<u64> = Vec::with_capacity(dictionary.len());
    dictionary.iter().rev().for_each(|matches| {
        if bases == 0 {
            bases = matches.0;
            let mut arr: [u8; 8] = [0; 8];
            arr[0..4].copy_from_slice(&(matches.1.start as u32).to_ne_bytes());
            eprintln!("{:?}", arr[0..4].to_vec());
            arr[4..5].copy_from_slice(&(matches.0 as u8).to_ne_bytes());
            arr[5..6].copy_from_slice(&(first as u8).to_ne_bytes());

            u64_encoding.push(u64::from_ne_bytes(arr));
            if first { first = false };
        }
        bases -= 1;
    });

    u64_encoding.shrink_to_fit();
    u64_encoding
}
