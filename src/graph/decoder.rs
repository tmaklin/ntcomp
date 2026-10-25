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
use super::header::FileHeader;
use crate::codecs::bits::BitReader;
use super::choices::Choices;
use super::snps::SnpEdits;
use crate::nruns::NRuns;
use crate::codecs::lzma::{self, LzmaData};
use super::ChainCsr;
use crate::walk::WalkIndex;

use crate::index::Index;

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

/// Reads `n_bytes` from `input`. The buffer grows with the bytes read, so a
/// corrupt size fails at the end of the input instead of allocating it.
fn read_bytes<R: Read>(
    input: &mut R,
    n_bytes: u64,
) -> Result<Vec<u8>, E> {
    let mut bytes = Vec::new();
    input.take(n_bytes).read_to_end(&mut bytes)?;
    if bytes.len() as u64 != n_bytes {
        return Err(Box::new(ArchiveError(format!("a section of {} bytes is truncated to {}", n_bytes, bytes.len()))))
    }
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

/// The sections of an archive ([GraphDecoder::decode_from_read]); the
/// optional ones are None if the archive does not have them.
pub struct DecodedSections {
    pub graph: ChainCsr,
    pub choices: Choices,
    /// Hash of every contig.
    pub hashes: Option<Vec<super::ContigHash>>,
    pub contig_names: Option<Vec<Vec<u8>>>,
    /// Input files: name and colours.
    pub file_ranges: Option<Vec<(Vec<u8>, Range<u32>)>>,
    /// Per contig: whether it is stored reverse complemented.
    pub reversed: Option<Vec<bool>>,
    pub snp_edits: Option<SnpEdits>,
    /// The runs of `N` cut out of the contigs.
    pub n_runs: Option<NRuns>,
}

/// Decodes the output from [GraphEncoder](super::encoder::GraphEncoder).
pub struct GraphDecoder<'a, R: Read> {
    reader: &'a mut R,
    file_header: FileHeader,
}

impl<'a, R: Read> GraphDecoder<'a, R> {
    /// Reads the file header from `reader`, see [FileHeader::read_from].
    pub fn new(
        reader: &'a mut R,
    ) -> Result<Self, E> {
        let file_header = FileHeader::read_from(reader)?;

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
    /// Fails for an archive written with an index, see
    /// [GraphDecoder::decode_from_read_with_index].
    pub fn decode_from_read(
        self,
    ) -> Result<DecodedSections, E> {
        self.decode_sections(None)
    }

    /// Same as [GraphDecoder::decode_from_read] for an archive written with
    /// [GraphEncoder::with_index](super::encoder::GraphEncoder::with_index)
    /// and `index`: an error before any section is read if the header names
    /// another index ([FileHeader::check_index]). (To decode without the LCS
    /// array, see [GraphDecoder::decode_from_read_with].)
    pub fn decode_from_read_with_index(
        self,
        index: &Index,
    ) -> Result<DecodedSections, E> {
        self.file_header.check_index(index.id())?;
        self.decode_sections(Some(index.walk_index()))
    }

    /// Same as [GraphDecoder::decode_from_read_with_index] for the index as a
    /// [WalkIndex], whose LCS array is optional (the values are then
    /// computed from the SBWT, slower). The SBWT, and the LCS array if
    /// given, are checked against the header first; their fingerprints are
    /// computed here (an [Index] keeps its own).
    pub fn decode_from_read_with(
        self,
        index: WalkIndex,
    ) -> Result<DecodedSections, E> {
        let sbwt = index.sbwt;
        self.file_header.check_index_parts(
            sbwt.k().try_into()?, sbwt.n_kmers() as u64,
            crate::index::sbwt_fingerprint(sbwt), index.lcs.map(crate::index::lcs_fingerprint))?;
        self.decode_sections(Some(index))
    }

    /// Reads and decodes every section with `index` (needed if the archive
    /// was written with one); _k_ comes from the header.
    fn decode_sections(
        self,
        index: Option<WalkIndex>,
    ) -> Result<DecodedSections, E> {
        let header = self.file_header;
        let reader = self.reader;
        if header.index.is_some() && index.is_none() {
            return Err(Box::new(ArchiveError("the archive was written with an index; decode it with the index".to_string())))
        }
        let k = header.index.map(|id| id.k);
        let n_contigs = header.n_contigs;

        let file_ranges = read_file_ranges(reader, header.file_range_bytes)?;
        if file_ranges.as_ref().is_some_and(|r| r.last().map_or(0, |(_, range)| range.end) != n_contigs) {
            return Err(Box::new(ArchiveError(format!("file ranges do not cover {} contigs", n_contigs))))
        }
        let contig_names = read_contig_names(reader, header.contig_name_bytes)?;
        if header.graph_bytes == 0 {
            return Err(Box::new(ArchiveError("the graph section is missing".to_string())))
        }
        let graph = ChainCsr::from_bytes(&read_bytes(reader, header.graph_bytes)?, index, k)?;
        let hashes: Option<Vec<super::ContigHash>> = read_section(reader, header.hash_bytes)?;
        if let Some(hashes) = &hashes && hashes.len() != n_contigs as usize {
            return Err(Box::new(ArchiveError(format!("{} hashes for {} contigs", hashes.len(), n_contigs))))
        }
        let reversed = read_reversed(reader, header.strand_bytes, n_contigs as usize)?;
        let snp_edits: Option<SnpEdits> = read_section(reader, header.snp_bytes)?;
        if snp_edits.as_ref().is_some_and(|edits| edits.n_colors() != n_contigs as usize) {
            return Err(Box::new(ArchiveError(format!("SNP section does not have {} colours", n_contigs))))
        }
        let n_runs = match header.n_run_bytes {
            0 => None,
            n_bytes => Some(NRuns::from_bytes(&read_bytes(reader, n_bytes)?, n_contigs as usize)?),
        };
        if header.choice_bytes == 0 {
            return Err(Box::new(ArchiveError("the choice section is missing".to_string())))
        }
        let choices = Choices::from_stored(&read_bytes(reader, header.choice_bytes)?, &graph, n_contigs as usize)?;

        Ok(DecodedSections { graph, choices, hashes, contig_names, file_ranges, reversed, snp_edits, n_runs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::encoder::GraphEncoder;
    use super::super::{GraphBuilder, GraphWalk};
    use super::super::header::FILE_HEADER_BYTES;
    use super::super::snps::Snp;
    use crate::testing::{dictionary, random_sequence, related_index};

    fn graph() -> (ChainCsr, Choices, Vec<Vec<(usize, Range<usize>)>>) {
        let dicts: Vec<Vec<(usize, Range<usize>)>> = vec![
            vec![(5, 10..11), (1, 50..51), (2, 60..61), (5, 10..11), (3, 20..21)],
            vec![(5, 10..11), (3, 20..21), (2, 40..41)],
        ];
        let mut builder = GraphBuilder::new();
        for d in &dicts {
            builder.add_dictionary(&dictionary(d));
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
        let runs = NRuns { runs: vec![vec![crate::nruns::NRun { pos: 2, len: 250 }], vec![]] };
        let bytes = encode(GraphEncoder::new(graph, choices)
            .with_contig_names(names.clone())
            .with_hashes(hashes.clone())
            .with_file_ranges(ranges.clone())
            .with_reversed(vec![false, true])
            .with_snp_edits(snps.clone())
            .with_n_runs(runs.clone()));

        let mut input = bytes.as_slice();
        let decoder = GraphDecoder::new(&mut input).unwrap();
        assert_eq!(decoder.header().n_contigs, 2);
        let DecodedSections { graph, choices, hashes: got_hashes, contig_names: got_names, file_ranges: got_ranges, reversed: got_reversed, snp_edits: got_snps, n_runs: got_runs } = decoder.decode_from_read().unwrap();
        assert_eq!(got_runs, Some(runs));
        assert_eq!(got_reversed, Some(vec![false, true]));
        assert_eq!(got_snps, Some(snps));
        assert_eq!(got_hashes, Some(hashes));
        assert_eq!(got_names, Some(names));
        assert_eq!(got_ranges, Some(ranges));
        for (color, d) in dicts.iter().enumerate() {
            let decoded = graph.walk_path(&choices, color as u32).unwrap();
            assert_eq!(decoded.sequence_len, d.iter().map(|x| x.0).sum::<usize>());
        }
    }

    #[test]
    fn absent_sections_are_none() {
        let (graph, choices, _) = graph();
        let bytes = encode(GraphEncoder::new(graph, choices));
        let mut input = bytes.as_slice();
        let sections = GraphDecoder::new(&mut input).unwrap().decode_from_read().unwrap();
        assert_eq!(sections.choices.streams.len(), 2);
        assert!(sections.hashes.is_none() && sections.contig_names.is_none() && sections.file_ranges.is_none() && sections.reversed.is_none() && sections.snp_edits.is_none());
    }

    #[test]
    fn strands_are_only_stored_if_a_contig_is_reversed() {
        let (stored, choices, _) = graph();
        let bytes = encode(GraphEncoder::new(stored, choices).with_reversed(vec![false, false]));
        let mut input = bytes.as_slice();
        let decoder = GraphDecoder::new(&mut input).unwrap();
        assert_eq!(decoder.header().strand_bytes, 0);
        assert!(decoder.decode_from_read().unwrap().reversed.is_none());

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
        use crate::graph::GraphWalk;
        let (seqs, index) = related_index(15);
        let mut builder = crate::graph::GraphBuilder::new();
        for seq in &seqs {
            builder.add_dictionary(&crate::nlz::parse(seq, &index).unwrap());
        }
        let (stored, choices) = builder.finish();
        let walks: Vec<_> = (0..seqs.len() as u32).map(|c| stored.walk_path(&choices, c).unwrap()).collect();
        let mut with_lcs = Vec::new();
        crate::graph::encoder::GraphEncoder::new(stored, choices).with_index(&index).encode_to_write(&mut with_lcs).unwrap();

        for array in [Some(index.lcs()), None] {
            let walk = crate::walk::WalkIndex { sbwt: index.sbwt(), lcs: array };
            let mut input = with_lcs.as_slice();
            let DecodedSections { graph, choices, .. } = crate::graph::decoder::GraphDecoder::new(&mut input).unwrap().decode_from_read_with(walk).unwrap();
            for (color, walk) in walks.iter().enumerate() {
                assert_eq!(&graph.walk_path(&choices, color as u32).unwrap(), walk);
            }
        }
        let mut input = with_lcs.as_slice();
        assert!(crate::graph::decoder::GraphDecoder::new(&mut input).unwrap().decode_from_read().is_err());
    }

    #[test]
    fn archives_name_their_index() {
        let (seqs, index) = related_index(15);
        let mut builder = GraphBuilder::new();
        for seq in &seqs {
            builder.add_dictionary(&crate::nlz::parse(seq, &index).unwrap());
        }
        let (graph, choices) = builder.finish();
        let mut bytes = Vec::new();
        GraphEncoder::new(graph, choices).with_index(&index).encode_to_write(&mut bytes).unwrap();

        let mut input = bytes.as_slice();
        let decoder = GraphDecoder::new(&mut input).unwrap();
        assert_eq!(decoder.header().index, Some(index.id()));
        assert!(decoder.decode_from_read_with_index(&index).is_ok());

        // Same k, other k-mers: refused before any section is read.
        let other = crate::index::Index::build(&[random_sequence(2000, 99)], &crate::index::BuildOpts { k: 15, ..Default::default() });
        assert_ne!(other.id().sbwt_fingerprint, index.id().sbwt_fingerprint);
        let err = GraphDecoder::new(&mut bytes.as_slice()).unwrap().decode_from_read_with_index(&other).err().unwrap();
        assert!(err.to_string().contains("different index"), "{}", err);
        // Also as a WalkIndex, with or without the LCS array.
        for lcs in [Some(other.lcs()), None] {
            let walk = crate::walk::WalkIndex { sbwt: other.sbwt(), lcs };
            assert!(GraphDecoder::new(&mut bytes.as_slice()).unwrap().decode_from_read_with(walk).is_err());
        }
        let walk = crate::walk::WalkIndex { sbwt: index.sbwt(), lcs: None };
        assert!(GraphDecoder::new(&mut bytes.as_slice()).unwrap().decode_from_read_with(walk).is_ok());
    }

    #[test]
    fn corrupt_section_sizes_fail_without_allocating() {
        let (graph, choices, _) = graph();
        let mut bytes = Vec::new();
        GraphEncoder::new(graph, choices).encode_to_write(&mut bytes).unwrap();
        // The graph section (third of eight sizes) claims 2^40 bytes.
        let at = FILE_HEADER_BYTES - 6 * 8;
        bytes[at..at + 8].copy_from_slice(&(1_u64 << 40).to_le_bytes());
        let err = GraphDecoder::new(&mut bytes.as_slice()).unwrap().decode_from_read().err().unwrap();
        assert!(err.to_string().contains("truncated"), "{}", err);
    }
}
