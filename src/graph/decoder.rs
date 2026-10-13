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
//! Decoder for reading NLZ-OG bytes.
//!
//! Reads the output from [GraphEncoder](super::encoder::GraphEncoder).
//!
use crate::{
    FileHeader,
    FILE_HEADER_BYTES,
    decode_file_header,
};
use super::choices::Choices;
use super::StoredGraph;

use core::ops::Range;

use std::io::Read;

use blake3::Hash;
use lzma_rust2::LzmaReader;

type E = Box<dyn std::error::Error>;


/// Helper for deserializing postcard data
fn read_section<T: serde::de::DeserializeOwned, R: Read>(
    input: &mut R,
    n_bytes: u64,
) -> Result<T, E> {
    let mut bytes = vec![0_u8; n_bytes as usize];
    input.read_exact(&mut bytes)?;
    let ret = postcard::from_bytes(&bytes)?;
    Ok(ret)
}

/// Decodes the output from [GraphEncoder](super::encoder::GraphEncoder).
pub struct GraphDecoder<'a, R: Read> {
    reader: &'a mut R,
    file_header: FileHeader,

    /// Hashes
    pub hashes: Vec<Hash>,

    /// Contig names
    pub contig_names: Vec<Vec<u8>>,

    /// File ranges
    pub file_ranges: Vec<(Vec<u8>, Range<u32>)>,

    /// Path lens
    pub path_lens: Vec<u32>,
}

impl<'a, R: Read> GraphDecoder<'a, R> {
    pub fn new(
        reader: &'a mut R,
    ) -> Result<Self, E> {
        let mut header_bytes = [0_u8; FILE_HEADER_BYTES];
        reader.read_exact(&mut header_bytes).unwrap();
        let file_header = decode_file_header(&header_bytes)?;

        Ok(GraphDecoder {
            reader,
            file_header,
            hashes: Vec::new(),
            contig_names: Vec::new(),
            file_ranges: Vec::new(),
            path_lens: Vec::new(),
        })
    }

    pub fn header(
        &self,
    ) -> FileHeader {
        self.file_header.clone()
    }

    pub fn decode_from_read(
        &mut self,
    ) -> Result<(StoredGraph, Choices), E> {
        self.file_ranges = read_section(&mut self.reader, self.file_header.file_range_bytes as u64)?;
        let mut lzma_bytes: Vec<u8> = vec![0; self.file_header.contig_name_bytes as usize];
        self.reader.read_exact(&mut lzma_bytes)?;
        self.contig_names = self.decode_contig_names(&lzma_bytes)?;
        let graph: StoredGraph = read_section(&mut self.reader, self.file_header.graph_bytes)?;
        self.hashes = read_section(&mut self.reader, self.file_header.hash_bytes as u64)?;
        self.path_lens = read_section(&mut self.reader, self.file_header.path_len_bytes as u64)?;
        let choices: Choices = read_section(&mut self.reader, self.file_header.choice_bytes)?;

        Ok((graph, choices))
    }

    fn decode_contig_names(
        &mut self,
        compressed: &[u8],
    ) -> Result<Vec<Vec<u8>>, E> {
        let mut reader =  LzmaReader::new_mem_limit(compressed, u32::MAX, None).unwrap();
        let mut buf = [0; 1024];
        let mut out = Vec::new();
        loop {
            let n = reader.read(&mut buf).unwrap();
            if n == 0 {
                break;
            }
            out.extend_from_slice(&buf[..n]);
        }

        let n_bytes = out.len();
        let contig_names: Vec<Vec<u8>> = read_section(&mut std::io::Cursor::new(out), n_bytes as u64)?;
        Ok(contig_names)
    }
}
