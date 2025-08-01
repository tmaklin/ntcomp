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
use dsi_bitstream::codes::MinimalBinaryRead;
use flate2::write::GzDecoder;

use crate::BlockHeader;
use crate::encode::Codec;

type E = Box<dyn std::error::Error>;

#[derive(Debug, Clone)]
struct DecodeError {
    message: String,
}

impl std::fmt::Display for DecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "ntx decoder error: {}", self.message)
    }
}

impl std::error::Error for DecodeError {}

fn inflate_bytes(
    deflated: &[u8],
) -> Result<Vec<u8>, E> {
    let mut inflated: Vec<u8> = Vec::new();
    let mut decoder = GzDecoder::new(&mut inflated);
    decoder.write_all(deflated)?;
    decoder.finish()?;
    Ok(inflated)
}

fn rice_decode(
    encoded: &[u64],
    n_records: usize,
    param: usize,
) -> Result<Vec<u64>, E> {
    let mut reader = BufBitReader::<BE, _>::new(MemWordReader::new(encoded));

    let mut decoded: Vec<u64> = vec![0; n_records];
    for i in 0..n_records {
        decoded[i] = reader.read_rice(param)?;
    }
    Ok(decoded)
}

fn minimal_binary_decode(
    encoded: &[u64],
    n_records: usize,
    param: u64,
) -> Result<Vec<u64>, E> {
    let mut reader = BufBitReader::<BE, _>::new(MemWordReader::new(encoded));

    let mut decoded: Vec<u64> = vec![0; n_records];
    for i in 0..n_records {
        decoded[i] = reader.read_minimal_binary(param)?;
    }
    Ok(decoded)
}

pub fn decompress_block(
    block: &[u8],
    header: &BlockHeader,
    codec: Codec,
) -> Result<Vec<u64>, E> {
    let bytes = inflate_bytes(block)?;
    let encoded: Vec<u64> = bytes.chunks(8).map(|x| {
        let arr: [u8; 8] = [x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]];
        u64::from_ne_bytes(arr)
    }).collect();

    if encoded.len() != header.encoded_size as usize {
        return Err(Box::new(DecodeError{ message: "encoded.len() != header.encoded_size".to_owned() }));
    }

    let res = match codec {
        Codec::Rice => rice_decode(&encoded, header.num_u64 as usize, header.rice_param as usize)?,
        Codec::MinimalBinary => minimal_binary_decode(&encoded, header.num_u64 as usize, header.rice_param as u64)?,
    };

    Ok(res)
}

pub fn zip_block_contents(
    colex_ranks: &[u64],
    match_lengths: &[u64],
    flags: &[u64],
    bitnuc_codings: &[u64],
) -> Result<Vec<u64>, E> {
    if colex_ranks.len() != match_lengths.len() {
        return Err(Box::new(DecodeError{ message: "1".to_owned() }));
    }
    if flags.len() - colex_ranks.len() != bitnuc_codings.len() {
        return Err(Box::new(DecodeError{ message: "2".to_owned() }));
    }

    let mut i = 0;
    let mut j = 0;
    let res = flags.iter().map(|x| {
        let flag_bytes = x.to_ne_bytes()[0];
        let mut arr: [u8; 8] = [0; 8];
        let flag = u8::from_ne_bytes([flag_bytes]);
        if flag & 0b00000010 == 0b00000000 {
            let colex_bytes: Vec<u8> = colex_ranks[i].to_ne_bytes()[0..4].to_vec();
            let matchlen_bytes: Vec<u8> = match_lengths[i].to_ne_bytes()[0..7].to_vec();
            arr[0..4].copy_from_slice(&colex_bytes[0..4]);
            arr[4..7].copy_from_slice(&matchlen_bytes[0..3]);
            i += 1;
        } else {
            let bitnuc_bytes: Vec<u8> = bitnuc_codings[j].to_ne_bytes()[0..7].to_vec();
            arr[0..7].copy_from_slice(&bitnuc_bytes[0..7]);
            j += 1;
        }
        arr[7..8].copy_from_slice(&[flag_bytes]);
        u64::from_ne_bytes(arr)
    }).collect::<Vec<u64>>();

    Ok(res)
}

pub fn decode_dictionary(
    encoded: &[u64],
) -> Vec<Vec<(u32, u32, bool)>> {
    let mut dictionary: Vec<(u32, u32, bool)> = Vec::new();
    let mut temp: Vec<(u32, u32, bool)> = Vec::new();
    encoded.iter().rev().for_each(|x| {
        let arr: [u8; 8] = x.to_ne_bytes();
        let mut arr1: [u8; 4] = [0; 4];
        let mut arr2: [u8; 4] = [0; 4];
        let mut arr3: [u8; 1] = [0; 1];
        arr1.copy_from_slice(&arr[0..4]);
        arr2[0..3].copy_from_slice(&arr[4..7]);
        arr3.copy_from_slice(&arr[7..8]);
        let colex_rank = u32::from_ne_bytes(arr1);
        let ms = u32::from_ne_bytes(arr2);
        temp.push((ms, colex_rank, u8::from_ne_bytes(arr3) == 1));
        if u8::from_ne_bytes(arr3) == 1 {
            dictionary.extend(temp.iter().rev());
            temp.clear();
        }
    });

    let mut dictionaries: Vec<Vec<(u32, u32, bool)>> = Vec::new();

    dictionary.iter().for_each(|record| {
        let n = dictionaries.len();
        if record.2 {
            dictionaries.push(vec![*record]);
        } else {
            dictionaries[n - 1].push(*record);
        }
    });
    dictionaries
}
