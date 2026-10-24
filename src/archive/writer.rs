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
//! Writing an archive.
//!
use std::io::Write;
use std::ops::Range;

use rayon::prelude::*;

use sbwt::LcsArray;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

use crate::fastx::Contig;
use crate::graph::encoder::GraphEncoder;
use crate::graph::{contig_hash, ContigHash, GraphBuilder};

type E = Box<dyn std::error::Error>;

/// An input file: its name and its contigs.
pub type InputFile = (Vec<u8>, Vec<Contig>);

/// Builds an archive from input files added in order. Each contig gets the
/// next colour; the archive keeps its hash, its name and its file.
///
/// With [GraphBuilder::edit_snps], every contig's strand and edits are
/// chosen in a first pass ([ArchiveWriter::call_edits]) over all the files
/// before they are added ([ArchiveWriter::add_files]) in the same order.
pub struct ArchiveWriter<'a> {
    builder: GraphBuilder,
    sbwt: &'a SbwtIndexVariant,
    lcs: &'a LcsArray,
    hashes: Vec<ContigHash>,
    contig_names: Vec<Vec<u8>>,
    file_ranges: Vec<(Vec<u8>, Range<u32>)>,
}

impl<'a> ArchiveWriter<'a> {
    /// A writer that adds the contigs to `builder` (with its options) and
    /// encodes them against the index `sbwt` with its LCS array `lcs`.
    pub fn new(
        builder: GraphBuilder,
        sbwt: &'a SbwtIndexVariant,
        lcs: &'a LcsArray,
    ) -> Self {
        ArchiveWriter { builder, sbwt, lcs, hashes: Vec::new(), contig_names: Vec::new(), file_ranges: Vec::new() }
    }

    /// First pass of editing: chooses the strands and calls the edits of the
    /// contigs of `files`, in parallel. Returns the colours they will get.
    pub fn call_edits(
        &mut self,
        files: &[InputFile],
    ) -> Vec<u32> {
        let sequences: Vec<&[u8]> = files.iter().flat_map(|(_, contigs)| contigs).map(|(_, sequence)| sequence.as_slice()).collect();
        self.builder.call_edits_batch(&sequences, self.sbwt, self.lcs)
    }

    /// Adds the contigs of `files`. Returns the colours they got.
    pub fn add_files(
        &mut self,
        files: &[InputFile],
    ) -> Result<Vec<u32>, E> {
        let sequences: Vec<&[u8]> = files.iter().flat_map(|(_, contigs)| contigs).map(|(_, sequence)| sequence.as_slice()).collect();
        self.hashes.par_extend(sequences.par_iter().map(|sequence| contig_hash(sequence)));
        let added = self.builder.add_contigs(&sequences, self.sbwt, self.lcs)?;
        let mut color = self.contig_names.len() as u32;
        debug_assert_eq!(added.first().map_or(color, |&c| c), color);
        for (file_name, contigs) in files {
            let start = color;
            self.contig_names.extend(contigs.iter().map(|(name, _)| name.clone()));
            color += contigs.len() as u32;
            self.file_ranges.push((file_name.clone(), start..color));
        }
        Ok(added)
    }

    /// Per colour added so far: whether the contig is stored reverse
    /// complemented.
    pub fn reversed(&self) -> &[bool] {
        self.builder.reversed()
    }

    /// Writes the archive to `out`. Returns a summary of the section sizes.
    pub fn write<W: Write>(
        self,
        out: &mut W,
    ) -> Result<String, E> {
        let reversed = self.builder.reversed().to_vec();
        let snp_edits = self.builder.snp_edits();
        let (graph, choices) = self.builder.finish();
        let mut encoder = GraphEncoder::new(graph, choices)
            .with_hashes(self.hashes)
            .with_contig_names(self.contig_names)
            .with_file_ranges(self.file_ranges)
            .with_reversed(reversed)
            .with_index(self.sbwt, self.lcs);
        if let Some(snp_edits) = snp_edits {
            encoder = encoder.with_snp_edits(snp_edits);
        }
        encoder.encode_to_write(out)
    }
}
