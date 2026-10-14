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
use super::bits::BitReader;
use super::choices::Choices;
use super::snps::SnpEdits;
use super::lzma_wrapper::{self, LzmaData};
use super::StoredGraph;

use core::ops::Range;

use std::io::Read;

use blake3::Hash;

type E = Box<dyn std::error::Error>;

#[derive(Debug)]
pub struct ArchiveError(pub String);

impl std::fmt::Display for ArchiveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "invalid archive: {}", self.0)
    }
}

impl std::error::Error for ArchiveError {}

/// Reads `n_bytes` from `input`.
fn read_bytes<R: Read>(
    input: &mut R,
    n_bytes: u64,
) -> Result<Vec<u8>, E> {
    let mut bytes = vec![0_u8; n_bytes.try_into()?];
    input.read_exact(&mut bytes)?;
    Ok(bytes)
}

/// Deserializes a section of `n_bytes` bytes, None if the section is empty (absent).
fn read_section<T: serde::de::DeserializeOwned, R: Read>(
    input: &mut R,
    n_bytes: u64,
) -> Result<Option<T>, E> {
    if n_bytes == 0 {
        return Ok(None)
    }
    Ok(Some(postcard::from_bytes(&read_bytes(input, n_bytes)?)?))
}

/// Deserializes a section that must be present.
fn read_required_section<T: serde::de::DeserializeOwned, R: Read>(
    input: &mut R,
    n_bytes: u64,
    name: &str,
) -> Result<T, E> {
    read_section(input, n_bytes)?.ok_or_else(|| Box::new(ArchiveError(format!("the {} section is missing", name))) as E)
}

/// Contig names, LZMA compressed; None if the section is empty (absent).
fn read_contig_names<R: Read>(
    input: &mut R,
    n_bytes: u64,
) -> Result<Option<Vec<Vec<u8>>>, E> {
    let Some(lzma) = read_section::<LzmaData, R>(input, n_bytes)? else {
        return Ok(None)
    };
    Ok(Some(postcard::from_bytes(&lzma_wrapper::decompress(&lzma)?)?))
}

/// One bit per colour, set if the contig is reverse complemented; None if
/// the section is empty (no contig is reversed).
fn read_reversed<R: Read>(
    input: &mut R,
    n_bytes: u64,
    n_colors: usize,
) -> Result<Option<Vec<bool>>, E> {
    let Some(bytes) = read_section::<Vec<u8>, R>(input, n_bytes)? else {
        return Ok(None)
    };
    if bytes.len() != n_colors.div_ceil(8) {
        return Err(Box::new(ArchiveError(format!("strand section has {} bytes for {} colours", bytes.len(), n_colors))))
    }
    let mut bits = BitReader::new(&bytes);
    Ok(Some((0..n_colors).map(|_| bits.read_bits(1) == Some(1)).collect()))
}

/// Sections of an archive, in the order returned by [GraphDecoder::decode_from_read].
pub type DecodedSections = (
    StoredGraph,
    Choices,
    Option<Vec<Hash>>,
    Option<Vec<Vec<u8>>>,
    Option<Vec<(Vec<u8>, Range<u32>)>>,
    Option<Vec<bool>>,
    Option<SnpEdits>,
);

/// Decodes the output from [GraphEncoder](super::encoder::GraphEncoder).
pub struct GraphDecoder<'a, R: Read> {
    reader: &'a mut R,
    file_header: FileHeader,
}

impl<'a, R: Read> GraphDecoder<'a, R> {
    /// Reads the file header from `reader`.
    pub fn new(
        reader: &'a mut R,
    ) -> Result<Self, E> {
        let mut header_bytes = [0_u8; FILE_HEADER_BYTES];
        reader.read_exact(&mut header_bytes)?;
        let file_header = decode_file_header(&header_bytes)?;

        Ok(GraphDecoder {
            reader,
            file_header,
        })
    }

    pub fn header(
        &self,
    ) -> FileHeader {
        self.file_header.clone()
    }

    /// Reads and decodes every section, see [DecodedSections].
    pub fn decode_from_read(
        self,
    ) -> Result<DecodedSections, E> {
        let header = self.file_header;
        let reader = self.reader;

        let file_ranges = read_section(reader, header.file_range_bytes as u64)?;
        let contig_names = read_contig_names(reader, header.contig_name_bytes)?;
        let graph: StoredGraph = read_required_section(reader, header.graph_bytes, "graph")?;
        let hashes = read_section(reader, header.hash_bytes as u64)?;
        let reversed = read_reversed(reader, header.strand_bytes as u64, header.n_queries as usize)?;
        let snp_edits: Option<SnpEdits> = read_section(reader, header.snp_bytes as u64)?;
        if snp_edits.as_ref().is_some_and(|edits| edits.n_colors() != header.n_queries as usize) {
            return Err(Box::new(ArchiveError(format!("SNP section does not have {} colours", header.n_queries))))
        }
        let choices: Choices = read_required_section(reader, header.choice_bytes, "choice")?;

        Ok((graph, choices, hashes, contig_names, file_ranges, reversed, snp_edits))
    }
}

    }
}
