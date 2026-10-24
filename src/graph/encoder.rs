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
use super::header::{FileHeader, FILE_HEADER_BYTES};
use super::StoredGraph;
use super::choices::Choices;
use super::snps::SnpEdits;
use crate::codecs::bits::BitWriter;
use crate::codecs::lzma;

use core::ops::Range;

use std::io::Write;

use postcard::to_allocvec;
use sbwt::LcsArray;
use crate::walk::WalkIndex;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

type E = Box<dyn std::error::Error>;

/// Encodes the output from [GraphBuilder](super::GraphBuilder).
pub struct GraphEncoder<'a> {
    /// Stored graph
    graph: StoredGraph,

    /// Choice stream
    choices: Choices,

    /// Hashes of the contigs, see [ContigHash](super::ContigHash)
    hashes: Option<Vec<super::ContigHash>>,

    /// Contig names
    contig_names: Option<Vec<Vec<u8>>>,

    /// File ranges
    file_ranges: Option<Vec<(Vec<u8>, Range<u32>)>>,

    /// Per colour: true if the contig is stored reverse complemented
    reversed: Option<Vec<bool>>,

    /// Edited SNPs
    snp_edits: Option<SnpEdits>,

    /// The index with its LCS array: the suffix lengths are coded per LCS
    /// context and graph nodes are derived by walks in it, see
    /// [collapsed](super::collapsed).
    index: Option<WalkIndex<'a>>,

    /// _k_-mer size of the index, for the direct suffix length symbols
    k: Option<u32>,
}

/// Serializes an optional section, 0 bytes if absent.
///
/// This is temporary and will be changed.
fn optional_section<T: serde::Serialize>(
    section: &Option<T>,
) -> Result<Vec<u8>, E> {
    match section {
        Some(value) => Ok(to_allocvec(value)?),
        None => Ok(Vec::new()),
    }
}

impl<'a> GraphEncoder<'a> {
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
            reversed: None,
            snp_edits: None,
            index: None,
            k: None,
        }
    }

    /// Add the contigs' hashes ([contig_hash](super::contig_hash))
    pub fn with_hashes(
        mut self,
        hashes: Vec<super::ContigHash>,
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

    /// Add the strands of the contigs, see [GraphBuilder::reversed](super::GraphBuilder::reversed).
    pub fn with_reversed(
        mut self,
        reversed: Vec<bool>,
    ) -> Self {
        self.reversed = Some(reversed);
        self
    }

    /// Code the graph for the index the paths were encoded with: the
    /// suffix lengths per LCS context of `lcs`, the direct suffix length
    /// symbols for the _k_-mer size of `sbwt`, and the nodes derived by walks
    /// in `sbwt`, see
    /// [collapsed](super::collapsed). Decoding the archive then needs the
    /// index, see [GraphDecoder::decode_from_read_with_index](super::decoder::GraphDecoder::decode_from_read_with_index).
    pub fn with_index(
        mut self,
        sbwt: &'a SbwtIndexVariant,
        lcs: &'a LcsArray,
    ) -> Self {
        let SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt;
        self.index = Some(WalkIndex { sbwt, lcs: Some(lcs) });
        self.k = Some(sbwt.k().try_into().expect("k does not fit in u32"));
        self
    }

    /// Add the edited SNPs, see [GraphBuilder::snp_edits](super::GraphBuilder::snp_edits).
    pub fn with_snp_edits(
        mut self,
        snp_edits: SnpEdits,
    ) -> Self {
        self.snp_edits = Some(snp_edits);
        self
    }

    /// Edited SNPs, one stream per colour; not written (0 bytes) if absent.
    fn encode_snp_edits(
        &self,
    ) -> Result<Vec<u8>, E> {
        if let Some(edits) = &self.snp_edits
            && edits.n_colors() != self.choices.streams.len() {
            return Err(format!("SNPs of {} colours for {} colours", edits.n_colors(), self.choices.streams.len()).into())
        }
        optional_section(&self.snp_edits)
    }

    /// One bit per colour, set if the contig is reverse complemented. Not
    /// written (0 bytes) if no contig is, which the decoder reads as all forward.
    fn encode_reversed(
        &self,
    ) -> Result<Vec<u8>, E> {
        match &self.reversed {
            Some(reversed) if reversed.iter().any(|&r| r) => {
                if reversed.len() != self.choices.streams.len() {
                    return Err(format!("{} strands for {} colours", reversed.len(), self.choices.streams.len()).into())
                }
                let mut bits = BitWriter::new();
                for &r in reversed {
                    bits.write_bits(r as u64, 1);
                }
                Ok(to_allocvec(&bits.finish())?)
            },
            _ => Ok(Vec::new()),
        }
    }

    /// File ranges as stored: the file names and the number of contigs in each
    /// file (the ranges are contiguous from colour 0), LZMA compressed.
    fn encode_file_ranges(
        &self,
    ) -> Result<Vec<u8>, E> {
        let Some(ranges) = &self.file_ranges else { return Ok(Vec::new()) };
        let mut next = 0;
        for (_, range) in ranges {
            if range.start != next || range.end < range.start {
                return Err("file ranges are not contiguous from colour 0".into())
            }
            next = range.end;
        }
        let names: Vec<&Vec<u8>> = ranges.iter().map(|(name, _)| name).collect();
        let counts: Vec<u32> = ranges.iter().map(|(_, range)| range.end - range.start).collect();
        Ok(to_allocvec(&lzma::compress(&to_allocvec(&(names, counts))?)?)?)
    }

    /// Contig names are compressed with LZMA.
    fn encode_contig_names(
        &self,
    ) -> Result<Vec<u8>, E> {
        match &self.contig_names {
            Some(contig_names) => Ok(to_allocvec(&lzma::compress(&to_allocvec(contig_names)?)?)?),
            None => Ok(Vec::new()),
        }
    }

    /// Writes the archive to `writer`.
    ///
    /// Returns a description of the graph and the size of every section for
    /// logging (subject to change).
    pub fn encode_to_write<W: Write>(
        self,
        writer: &mut W,
    ) -> Result<String, E> {
        let n_colors: u32 = self.choices.streams.len().try_into()?;

        let file_range_bytes = self.encode_file_ranges()?;
        let contig_name_bytes = self.encode_contig_names()?;
        let hash_bytes = optional_section(&self.hashes)?;
        let strand_bytes = self.encode_reversed()?;
        let snp_bytes = self.encode_snp_edits()?;
        let graph_bytes = self.graph.to_bytes(self.index, self.k)?;
        let choice_bytes = self.choices.to_stored(&self.graph)?;

        let header = FileHeader{
            nlz_header: [0_u8; 6],
            n_queries: n_colors,
            contig_name_bytes: contig_name_bytes.len() as u64,
            graph_bytes: graph_bytes.len() as u64,
            hash_bytes: hash_bytes.len().try_into()?,
            strand_bytes: strand_bytes.len().try_into()?,
            snp_bytes: snp_bytes.len().try_into()?,
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
        writer.write_all(&strand_bytes)?;
        writer.write_all(&snp_bytes)?;
        writer.write_all(&choice_bytes)?;

        let total = FILE_HEADER_BYTES + file_range_bytes.len() + contig_name_bytes.len() + graph_bytes.len()
            + hash_bytes.len() + strand_bytes.len() + snp_bytes.len() + choice_bytes.len();
        let n_reversed = self.reversed.as_ref().map_or(0, |r| r.iter().filter(|&&r| r).count());
        let n_snps = self.snp_edits.as_ref().map_or(0, |s| s.n_snps());
        Ok(format!(
            "Graph: {}, {} bytes; choices: {} choices, {} bytes; contig names {} bytes, hashes {} bytes, strands {} bytes ({} reversed), SNPs {} bytes ({} edited), file ranges {} bytes; total {} bytes",
            self.graph.summary(), graph_bytes.len(), self.choices.n_choices(), choice_bytes.len(),
            contig_name_bytes.len(), hash_bytes.len(), strand_bytes.len(), n_reversed, snp_bytes.len(), n_snps, file_range_bytes.len(), total,
        ))
    }
}
