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
//! The NLZ-OG archive header: the sizes of the sections that follow it.
//!
use bincode::{Encode, Decode};
use bincode::decode_from_slice;

type E = Box<dyn std::error::Error>;

/// Header of an archive written by [GraphEncoder](super::encoder::GraphEncoder).
#[derive(Clone, Debug, Decode, Encode, PartialEq)]
pub struct FileHeader {
    /// Nlz header, consists of four ASCII bytes spelling "nrlz" and two bytes specifying the version.
    ///
    /// First four bytes can be used to check that a binary record is an .nlz record.
    /// Next two bytes can be used to check which version of nlz was used to generate this file.
    pub nlz_header: [u8; 6], // = [110, 114, 108, 122, ...];

    pub file_range_bytes: u32,

    pub hash_bytes: u32,

    /// Size of the section that marks the reverse complemented contigs.
    pub strand_bytes: u32,

    /// Size of the section with the edited SNPs.
    pub snp_bytes: u32,

    pub n_queries: u32,

    pub contig_name_bytes: u64,

    pub graph_bytes: u64,

    pub choice_bytes: u64,
}

pub const FILE_HEADER_BYTES: usize = 50;

/// Reads the header from its [FILE_HEADER_BYTES] bytes.
pub fn decode_file_header(
    header_bytes: &[u8],
) -> Result<FileHeader, E> {
    Ok(decode_from_slice(header_bytes, bincode::config::standard().with_fixed_int_encoding())?.0)
}

