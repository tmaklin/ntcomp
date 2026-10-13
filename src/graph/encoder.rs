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
//! Encoder for creating a valid NLZ-OG archive.
//!
//! Can be read back with [GraphDecoder](super::decoder::GraphDecoder).
use crate::FileHeader;
use crate::FILE_HEADER_BYTES;
use super::StoredGraph;
use super::choices::Choices;

use core::ops::Range;

use std::io::Write;

use blake3::Hash;
use lzma_rust2::{LzmaOptions, LzmaWriter};
use postcard::to_allocvec;

type E = Box<dyn std::error::Error>;

/// Encodes the output from [GraphBuilder](super::GraphBuilder).
pub struct GraphEncoder {
    /// Stored graph
    graph: StoredGraph,

    /// Choice stream
    choices: Choices,

    /// Hashes
    hashes: Option<Vec<Hash>>,

    /// Contig names
    contig_names: Option<Vec<Vec<u8>>>,

    /// File ranges
    file_ranges: Option<Vec<(Vec<u8>, Range<u32>)>>,

    /// Path lens
    path_lens: Option<Vec<u32>>,

}

impl GraphEncoder {
    pub fn new(
        graph: StoredGraph,
        choices: Choices,
    ) -> Self {
        GraphEncoder {
            graph,
            choices,
            hashes: None,
            contig_names: None,
            file_ranges: None,
            path_lens: None,
        }
    }

    /// Add hashes
    pub fn with_hashes(
        mut self,
        hashes: Vec<Hash>,
    ) -> Self {
        self.hashes = Some(hashes);
        self
    }

    /// Add contig_names
    pub fn with_contig_names(
        mut self,
        contig_names: Vec<Vec<u8>>,
    ) -> Self {
        self.contig_names = Some(contig_names);
        self
    }

    /// Add file ranges
    pub fn with_file_ranges(
        mut self,
        file_ranges: Vec<(Vec<u8>, Range<u32>)>,
    ) -> Self {
        self.file_ranges = Some(file_ranges);
        self
    }

    /// Add path lens
    pub fn with_path_lens(
        mut self,
        path_lens: Vec<u32>,
    ) -> Self {
        self.path_lens = Some(path_lens);
        self
    }

    pub fn summary(
        &self,
    ) -> Result<String, E> {
        let graph_bytes = to_allocvec(&self.graph)?;
        let mut ret = format!("Graph: {}, {} bytes", self.graph.summary(), graph_bytes.len());
        let choice_bytes = to_allocvec(&self.choices)?;
        ret += &format!("Choices: {} choices, {} bytes", self.choices.n_choices(), choice_bytes.len());
        Ok(ret)
    }

    fn encode_file_ranges(
        &mut self,
    ) -> Option<Result<Vec<u8>, E>> {
        let file_ranges = std::mem::take(self.file_ranges.as_mut()?);
        self.file_ranges = None;
        match to_allocvec(&file_ranges) {
            Ok(bytes) => Some(Ok(bytes)),
            Err(e) => Some(Err(Box::new(e))),
        }
    }

    fn encode_contig_names(
        &mut self,
    ) -> Option<Result<Vec<u8>, E>> {
        let contig_names = std::mem::take(self.contig_names.as_mut()?);
        self.contig_names = None;
        match to_allocvec(&contig_names) {
            Ok(bytes) => {
                let mut options = LzmaOptions::with_preset(9);
                options.lc = 0;
                options.lp = 0;
                options.pb = 0;
                options.dict_size = options.dict_size.min((bytes.len() as u32).max(4096));

                let mut writer = LzmaWriter::new_use_header(Vec::new(), &options, None).expect("creating LZMA writer failed");
                writer.write_all(&bytes).expect("LZMA compression failed");
                let lzma = writer.finish().expect("LZMA compression failed");

                Some(Ok(lzma))
            },
            Err(e) => Some(Err(Box::new(e))),
        }
    }

    fn encode_hashes(
        &mut self,
    ) -> Option<Result<Vec<u8>, E>> {
        let hashes = std::mem::take(self.hashes.as_mut()?);
        self.hashes = None;
        match to_allocvec(&hashes) {
            Ok(bytes) => Some(Ok(bytes)),
            Err(e) => Some(Err(Box::new(e))),
        }
    }

    fn encode_path_lens(
        &mut self,
    ) -> Option<Result<Vec<u8>, E>> {
        let path_lens = std::mem::take(self.path_lens.as_mut()?);
        self.path_lens = None;
        match to_allocvec(&path_lens) {
            Ok(bytes) => Some(Ok(bytes)),
            Err(e) => Some(Err(Box::new(e))),
        }
    }

    pub fn encode_to_write<W: Write>(
        &mut self,
        writer: &mut W,
    ) -> Result<(), E> {
        let n_colors: u32 = self.path_lens.as_ref().unwrap().len().try_into()?;

        let file_range_bytes = if let Some(bytes) = self.encode_file_ranges() {
            bytes?
        } else {
            Vec::new()
        };

        let contig_name_bytes = if let Some(bytes) = self.encode_contig_names() {
            bytes?
        } else {
            Vec::new()
        };

        let hash_bytes = if let Some(bytes) = self.encode_hashes() {
            bytes?
        } else {
            Vec::new()
        };

        let path_len_bytes = if let Some(bytes) = self.encode_path_lens() {
            bytes?
        } else {
            Vec::new()
        };

        let graph_bytes = to_allocvec(&self.graph)?;
        let choice_bytes = to_allocvec(&self.choices)?;

        let header = FileHeader{
            nlz_header: [0_u8; 6],
            path_len_bytes: path_len_bytes.len().try_into()?,
            n_queries: n_colors,
            contig_name_bytes: contig_name_bytes.len() as u64,
            graph_bytes: graph_bytes.len() as u64,
            hash_bytes: hash_bytes.len().try_into()?,
            file_range_bytes: file_range_bytes.len().try_into()?,
            choice_bytes: choice_bytes.len() as u64,
        };
        let nbytes = bincode::encode_into_std_write(
            &header,
            writer,
            bincode::config::standard().with_fixed_int_encoding(),
        )?;

        assert_eq!(nbytes, FILE_HEADER_BYTES);

        writer.write_all(&file_range_bytes)?;
        writer.write_all(&contig_name_bytes)?;
        writer.write_all(&graph_bytes)?;
        writer.write_all(&hash_bytes)?;
        writer.write_all(&path_len_bytes)?;
        writer.write_all(&choice_bytes)?;

        Ok(())
    }
}
