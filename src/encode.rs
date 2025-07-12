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

use bincode::encode_into_std_write;
use bincode::{Encode, Decode};
use flate2::write::GzEncoder;
use flate2::Compression;

use dsi_bitstream::impls::MemWordWriterVec;
use dsi_bitstream::impls::BufBitWriter;
use dsi_bitstream::traits::BE;
use dsi_bitstream::codes::RiceWrite;

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

pub fn encode_block(
    u64_encoding: &[u64],
    num_records: usize,
) -> Vec<u8> {
    // Best: Rice
    let inv_mean: f64 = ((u64_encoding.len() as f64).ln() - (u64_encoding.iter().sum::<u64>() as f64).ln()).exp();
    let param = dsi_bitstream::codes::rice::log2_b(inv_mean);

    let word_write = MemWordWriterVec::new(Vec::<u64>::new());
    let mut writer = BufBitWriter::<BE, _>::new(word_write);

    u64_encoding.iter().for_each(|n| { writer.write_rice(*n, param).unwrap(); } );
    let _ = writer.flush();
    let rice_encoded = writer.into_inner().unwrap().into_inner();

    let bytes = rice_encoded.iter().flat_map(|x| {
        x.to_ne_bytes()
    }).collect::<Vec<u8>>();

    let mut deflated: Vec<u8> = Vec::with_capacity(bytes.len());
    let mut encoder = GzEncoder::new(&mut deflated, Compression::default());
    encoder.write_all(&bytes).unwrap();
    encoder.finish().unwrap();

    let block_header = BlockHeader{ block_size: deflated.len() as u32,
                                    num_records: num_records as u32,
                                    num_u64: u64_encoding.len() as u32,
                                    encoded_size: rice_encoded.len() as u32,
                                    rice_param: param as u8,
                                    bitpacker_exponent: 8_u8,
                                    placeholder1: 0, placeholder2: 0, placeholder3: 0,
    };

    let mut block: Vec<u8> = Vec::with_capacity(32 + deflated.len()/8);

    let nbytes = encode_into_std_write(
        &block_header,
        &mut block,
        bincode::config::standard().with_fixed_int_encoding(),
    );
    assert_eq!(nbytes.unwrap(), 32);
    block.append(&mut deflated);
    block.shrink_to_fit();

    block
}
