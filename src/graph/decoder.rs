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
use super::header::{decode_file_header, FileHeader, FILE_HEADER_BYTES};
use crate::codecs::bits::BitReader;
use super::choices::Choices;
use super::snps::SnpEdits;
use crate::codecs::lzma::{self, LzmaData};
use super::StoredGraph;
use crate::walk::WalkIndex;

use sbwt::LcsArray;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

use core::ops::Range;

use std::io::Read;


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

/// Contig names, LZMA compressed; None if the section is empty (absent).
fn read_contig_names<R: Read>(
    input: &mut R,
    n_bytes: u64,
) -> Result<Option<Vec<Vec<u8>>>, E> {
    let Some(section) = read_section::<LzmaData, R>(input, n_bytes)? else {
        return Ok(None)
    };
    Ok(Some(postcard::from_bytes(&lzma::decompress(&section)?)?))
}

/// Reads the file ranges written by the encoder (file names and contig counts,
/// LZMA compressed); the ranges are contiguous from colour 0.
fn read_file_ranges<R: Read>(
    input: &mut R,
    n_bytes: u64,
) -> Result<Option<Vec<(Vec<u8>, Range<u32>)>>, E> {
    let Some(section) = read_section::<LzmaData, R>(input, n_bytes)? else {
        return Ok(None)
    };
    let (names, counts): (Vec<Vec<u8>>, Vec<u32>) = postcard::from_bytes(&lzma::decompress(&section)?)?;
    if names.len() != counts.len() {
        return Err(Box::new(ArchiveError("file ranges: names and counts differ".to_string())))
    }
    let mut start = 0_u32;
    let mut ranges = Vec::with_capacity(names.len());
    for (name, count) in names.into_iter().zip(counts) {
        let end = start.checked_add(count).ok_or_else(|| ArchiveError("file ranges: too many contigs".to_string()))?;
        ranges.push((name, start..end));
        start = end;
    }
    Ok(Some(ranges))
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
    Option<Vec<super::ContigHash>>,
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
    ///
    /// Fails for an archive written with the LCS array, see
    /// [GraphDecoder::decode_from_read_with].
    pub fn decode_from_read(
        self,
    ) -> Result<DecodedSections, E> {
        self.decode_from_read_with(None, None)
    }

    /// Same as [GraphDecoder::decode_from_read] for an archive written with
    /// [GraphEncoder::with_index](super::encoder::GraphEncoder::with_index):
    /// the _k_-mer size comes from `sbwt`, and the LCS values from `lcs` if
    /// given, else they are computed from `sbwt` (slower).
    pub fn decode_from_read_with_index(
        self,
        sbwt: &SbwtIndexVariant,
        lcs: Option<&LcsArray>,
    ) -> Result<DecodedSections, E> {
        let SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt;
        self.decode_from_read_with(Some(WalkIndex { sbwt, lcs }), Some(sbwt.k().try_into()?))
    }

    /// Same as [GraphDecoder::decode_from_read]; `index` is needed if the
    /// archive was written with it (LCS contexts and graph nodes derived by
    /// walks; its LCS array is optional), and `k` is the _k_-mer size it was
    /// written with (None works for _k_ ≤ 255).
    pub fn decode_from_read_with(
        self,
        index: Option<WalkIndex>,
        k: Option<u32>,
    ) -> Result<DecodedSections, E> {
        let header = self.file_header;
        let reader = self.reader;

        let file_ranges = read_file_ranges(reader, header.file_range_bytes as u64)?;
        if file_ranges.as_ref().is_some_and(|r| r.last().map_or(0, |(_, range)| range.end) != header.n_queries) {
            return Err(Box::new(ArchiveError(format!("file ranges do not cover {} contigs", header.n_queries))))
        }
        let contig_names = read_contig_names(reader, header.contig_name_bytes)?;
        if header.graph_bytes == 0 {
            return Err(Box::new(ArchiveError("the graph section is missing".to_string())))
        }
        let graph = StoredGraph::from_bytes(&read_bytes(reader, header.graph_bytes)?, index, k)?;
        let hashes = read_section(reader, header.hash_bytes as u64)?;
        let reversed = read_reversed(reader, header.strand_bytes as u64, header.n_queries as usize)?;
        let snp_edits: Option<SnpEdits> = read_section(reader, header.snp_bytes as u64)?;
        if snp_edits.as_ref().is_some_and(|edits| edits.n_colors() != header.n_queries as usize) {
            return Err(Box::new(ArchiveError(format!("SNP section does not have {} colours", header.n_queries))))
        }
        if header.choice_bytes == 0 {
            return Err(Box::new(ArchiveError("the choice section is missing".to_string())))
        }
        let choices = Choices::from_stored(&read_bytes(reader, header.choice_bytes)?, &graph, header.n_queries as usize)?;

        Ok((graph, choices, hashes, contig_names, file_ranges, reversed, snp_edits))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::encoder::GraphEncoder;
    use super::super::{GraphBuilder, GraphWalk};
    use super::super::snps::Snp;
    use crate::testing::related_index;

    fn graph() -> (StoredGraph, Choices, Vec<Vec<(usize, Range<usize>)>>) {
        let dicts: Vec<Vec<(usize, Range<usize>)>> = vec![
            vec![(5, 10..11), (1, 50..51), (2, 60..61), (5, 10..11), (3, 20..21)],
            vec![(5, 10..11), (3, 20..21), (2, 40..41)],
        ];
        let mut builder = GraphBuilder::new().collapse_chains(true);
        for d in &dicts {
            builder.add_path(d.iter().map(|x| x.0).sum(), d);
        }
        let (graph, choices) = builder.finish();
        (graph, choices, dicts)
    }

    fn encode(encoder: GraphEncoder) -> Vec<u8> {
        let mut bytes = Vec::new();
        encoder.encode_to_write(&mut bytes).unwrap();
        bytes
    }

    #[test]
    fn archive_roundtrip() {
        let (graph, choices, dicts) = graph();
        let names = vec![b"contig one".to_vec(), b"contig two".to_vec()];
        let hashes = vec![super::super::contig_hash(b"one"), super::super::contig_hash(b"two")];
        let ranges = vec![(b"file".to_vec(), 0..2)];
        let snps = SnpEdits::new(&[vec![], vec![Snp { pos: 3, base: b'A', index_base: b'G' }]]);
        let bytes = encode(GraphEncoder::new(graph, choices)
            .with_contig_names(names.clone())
            .with_hashes(hashes.clone())
            .with_file_ranges(ranges.clone())
            .with_reversed(vec![false, true])
            .with_snp_edits(snps.clone()));

        let mut input = bytes.as_slice();
        let decoder = GraphDecoder::new(&mut input).unwrap();
        assert_eq!(decoder.header().n_queries, 2);
        let (graph, choices, got_hashes, got_names, got_ranges, got_reversed, got_snps) = decoder.decode_from_read().unwrap();
        assert_eq!(got_reversed, Some(vec![false, true]));
        assert_eq!(got_snps, Some(snps));
        assert_eq!(got_hashes, Some(hashes));
        assert_eq!(got_names, Some(names));
        assert_eq!(got_ranges, Some(ranges));
        for (color, d) in dicts.iter().enumerate() {
            let (len, _) = graph.walk_path(&choices, color as u32).unwrap();
            assert_eq!(len as usize, d.iter().map(|x| x.0).sum::<usize>());
        }
    }

    #[test]
    fn absent_sections_are_none() {
        let (graph, choices, _) = graph();
        let bytes = encode(GraphEncoder::new(graph, choices));
        let mut input = bytes.as_slice();
        let (_, choices, hashes, names, ranges, reversed, snps) = GraphDecoder::new(&mut input).unwrap().decode_from_read().unwrap();
        assert_eq!(choices.streams.len(), 2);
        assert!(hashes.is_none() && names.is_none() && ranges.is_none() && reversed.is_none() && snps.is_none());
    }

    #[test]
    fn strands_are_only_stored_if_a_contig_is_reversed() {
        let (stored, choices, _) = graph();
        let bytes = encode(GraphEncoder::new(stored, choices).with_reversed(vec![false, false]));
        let mut input = bytes.as_slice();
        let decoder = GraphDecoder::new(&mut input).unwrap();
        assert_eq!(decoder.header().strand_bytes, 0);
        assert!(decoder.decode_from_read().unwrap().5.is_none());

        // One strand per colour is required.
        let (stored, choices, _) = graph();
        let mut bytes = Vec::new();
        assert!(GraphEncoder::new(stored, choices).with_reversed(vec![true]).encode_to_write(&mut bytes).is_err());
    }

    #[test]
    fn snps_need_one_stream_per_colour() {
        let (stored, choices, _) = graph();
        let mut bytes = Vec::new();
        let snps = SnpEdits::new(&[vec![Snp { pos: 0, base: b'C', index_base: b'T' }]]);
        assert!(GraphEncoder::new(stored, choices).with_snp_edits(snps).encode_to_write(&mut bytes).is_err());
    }

    #[test]
    fn truncated_archives_are_errors() {
        let (graph, choices, _) = graph();
        let bytes = encode(GraphEncoder::new(graph, choices).with_contig_names(vec![b"a".to_vec(), b"b".to_vec()]));
        for len in [0, 10, FILE_HEADER_BYTES, bytes.len() / 2, bytes.len() - 1] {
            let mut input = &bytes[..len];
            let decoded = GraphDecoder::new(&mut input).and_then(|decoder| decoder.decode_from_read());
            assert!(decoded.is_err(), "archive truncated to {} bytes", len);
        }
    }

    #[test]
    fn missing_required_section_is_an_error() {
        let (graph, choices, _) = graph();
        let mut bytes = encode(GraphEncoder::new(graph, choices));
        // Set the size of the choice section (last header field) to 0.
        bytes[FILE_HEADER_BYTES - 8..FILE_HEADER_BYTES].copy_from_slice(&0_u64.to_le_bytes());
        let mut input = bytes.as_slice();
        assert!(GraphDecoder::new(&mut input).unwrap().decode_from_read().is_err());
    }

    #[test]
    fn archives_written_with_the_index_decode_with_or_without_its_lcs_array() {
        use crate::graph::{StoredGraph, GraphWalk};
        let (seqs, variant, lcs) = related_index(15);
        let mut builder = crate::graph::GraphBuilder::new().collapse_chains(true);
        for seq in &seqs {
            builder.add_contig(seq, &variant, &lcs).unwrap();
        }
        let (stored, choices) = builder.finish();
        let walks: Vec<_> = (0..seqs.len() as u32).map(|c| stored.walk_path(&choices, c).unwrap()).collect();
        let mut with_lcs = Vec::new();
        crate::graph::encoder::GraphEncoder::new(stored, choices).with_index(&variant, &lcs).encode_to_write(&mut with_lcs).unwrap();

        let SbwtIndexVariant::SubsetMatrix(sbwt) = &variant;
        for array in [Some(&lcs), None] {
            let index = crate::walk::WalkIndex { sbwt, lcs: array };
            let mut input = with_lcs.as_slice();
            let (graph, choices, ..) = crate::graph::decoder::GraphDecoder::new(&mut input).unwrap().decode_from_read_with(Some(index), Some(15)).unwrap();
            for (color, walk) in walks.iter().enumerate() {
                assert_eq!(&graph.walk_path(&choices, color as u32).unwrap(), walk);
            }
            assert!(matches!(graph, StoredGraph::Chains(_)));
        }
        let mut input = with_lcs.as_slice();
        assert!(crate::graph::decoder::GraphDecoder::new(&mut input).unwrap().decode_from_read().is_err());
    }
}
