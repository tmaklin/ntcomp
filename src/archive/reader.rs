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
//! Reading an archive: any of its contigs, in any order.
//!
use std::collections::HashMap;
use std::io::Read;
use std::ops::Range;

use rayon::prelude::*;

use crate::graph::decoder::{ArchiveError, DecodedSections, GraphDecoder};
use crate::graph::snps::{self, SnpEdits};
use crate::graph::{contig_digest, file_hash, Choices, ChainCsr, ContigDigest, Csr, FileHash, GraphWalk};
use crate::index::Index;
use crate::nlz::PhraseCache;
use crate::nruns::{self, NRuns};

type E = Box<dyn std::error::Error>;

/// Contexts of SNPs in edited sequences, by colour and position, for putting
/// back repeated SNPs (see [snps]).
type SnpContexts = HashMap<(u32, u32), [u8; snps::CONTEXT]>;

/// An archive read with the index it was written with.
pub struct ArchiveReader<'a> {
    index: &'a Index,
    /// The phrases decoded so far.
    cache: PhraseCache,
    n_contigs: u32,
    graph: ChainCsr,
    choices: Choices,
    /// Per input file (one for all contigs without file ranges).
    pub(crate) hashes: Option<Vec<FileHash>>,
    contig_names: Option<Vec<Vec<u8>>>,
    file_ranges: Option<Vec<(Vec<u8>, Range<u32>)>>,
    reversed: Option<Vec<bool>>,
    snp_edits: Option<SnpEdits>,
    n_runs: Option<NRuns>,
}

/// A decoded contig.
pub struct DecodedContig {
    pub color: u32,
    pub sequence: Vec<u8>,
}

/// The check of an input file whose contigs have all been decoded
/// ([ContigDecoder::finished_files]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FileCheck {
    /// The file's index (0 without file ranges).
    pub file: usize,
    /// Whether its contigs match the stored hash.
    pub hash_match: bool,
}

impl<'a> ArchiveReader<'a> {
    /// Reads the archive from `input`; `index` is the index it was written
    /// with.
    pub fn open<R: Read>(
        input: &mut R,
        index: &'a Index,
    ) -> Result<Self, E> {
        let decoder = GraphDecoder::new(input)?;
        let n_contigs = decoder.header().n_contigs;
        let DecodedSections { graph, choices, hashes, contig_names, file_ranges, reversed, snp_edits, n_runs } = decoder.decode_from_read_with_index(index)?;
        if contig_names.as_ref().is_some_and(|names| names.len() != n_contigs as usize) {
            return Err(Box::new(ArchiveError(format!("{} contig names for {} contigs", contig_names.map_or(0, |n| n.len()), n_contigs))))
        }
        Ok(ArchiveReader { index, cache: PhraseCache::default(), n_contigs, graph, choices, hashes, contig_names, file_ranges, reversed, snp_edits, n_runs })
    }

    /// Number of contigs (colours).
    pub fn n_contigs(&self) -> u32 {
        self.n_contigs
    }

    /// Number of input files (0 if the archive does not store them).
    pub fn n_files(&self) -> usize {
        self.file_ranges.as_ref().map_or(0, |ranges| ranges.len())
    }

    /// Name of contig `color` (None if the archive has no names).
    pub fn contig_name(
        &self,
        color: u32,
    ) -> Option<&[u8]> {
        self.contig_names.as_ref().and_then(|names| names.get(color as usize)).map(|name| name.as_slice())
    }

    /// Whether the archive stores hashes to check decoded files by.
    pub fn has_hashes(&self) -> bool {
        self.hashes.is_some()
    }

    /// The input file of contig `color`: its index (0 without file ranges).
    pub fn file_of(
        &self,
        color: u32,
    ) -> usize {
        self.file_ranges.as_ref().map_or(0, |ranges| ranges.partition_point(|(_, range)| range.end <= color))
    }

    /// The contigs of input file `file` (all without file ranges).
    pub fn file_colors(
        &self,
        file: usize,
    ) -> Range<u32> {
        self.file_ranges.as_ref().map_or(0..self.n_contigs, |ranges| ranges.get(file).map_or(0..0, |(_, range)| range.clone()))
    }

    /// The colours of the contigs of the input file named `name`. An error if
    /// the archive has no file names or none matches.
    pub fn find_file(
        &self,
        name: &str,
    ) -> Result<Vec<u32>, E> {
        let ranges = self.file_ranges.as_ref().ok_or("the archive has no file names")?;
        let (_, range) = ranges.iter().find(|(file, _)| file.as_slice() == name.as_bytes())
            .ok_or_else(|| format!("file {} not found", name))?;
        Ok(range.clone().collect())
    }

    /// Name of the input file of contig `color` (None if the archive does
    /// not store the files).
    pub fn file_name(
        &self,
        color: u32,
    ) -> Option<&[u8]> {
        let ranges = self.file_ranges.as_ref()?;
        ranges.iter().find(|(_, range)| range.contains(&color)).map(|(name, _)| name.as_slice())
    }

    /// The colours of the contigs named `name`: the whole name or its first
    /// word. An error if the archive has no names or none matches.
    pub fn find_contig(
        &self,
        name: &str,
    ) -> Result<Vec<u32>, E> {
        let names = self.contig_names.as_ref().ok_or("the archive has no contig names")?;
        let name = name.as_bytes();
        let matches: Vec<u32> = names.iter().enumerate()
            .filter(|(_, x)| x.as_slice() == name || x.split(|c| *c == b' ').next() == Some(name))
            .map(|(i, _)| i as u32)
            .collect();
        if matches.is_empty() {
            return Err(format!("contig {} not found", String::from_utf8_lossy(name)).into())
        }
        Ok(matches)
    }

    /// A decoder for the contigs `colors`, in chunks of them in this order
    /// ([ContigDecoder::decode]). An error if a colour is out of range.
    pub fn decoder(
        &self,
        colors: &[u32],
    ) -> Result<ContigDecoder<'_, 'a>, E> {
        if let Some(bad) = colors.iter().find(|&&c| c >= self.n_contigs) {
            return Err(format!("color {} out of range ({} contigs)", bad, self.n_contigs).into())
        }
        // Repeated SNPs are located by contexts from the edited sequences of
        // earlier colours: the positions needed, by colour.
        let mut wanted: HashMap<u32, Vec<u32>> = HashMap::new();
        if let Some(edits) = &self.snp_edits {
            for &color in colors {
                for (source, pos) in edits.context_sources(color) {
                    wanted.entry(source).or_default().push(pos);
                }
            }
        }
        Ok(ContigDecoder { reader: self, wanted, contexts: SnpContexts::new(), pending: HashMap::new(), finished: Vec::new() })
    }

    /// Decodes the stored (edited) sequence of `color` from the graph,
    /// before the SNPs are put back.
    fn decode_edited(
        &self,
        color: u32,
    ) -> Result<Vec<u8>, String> {
        let dictionary = self.graph.walk_path(&self.choices, color).map_err(|e| e.to_string())?;
        dictionary.decode_cached(self.index, &self.cache).map_err(|e| e.to_string())
    }

    /// Puts the SNPs of `color` back into its edited sequence, then its runs
    /// of `N`, and restores its strand.
    fn finish_color(
        &self,
        color: u32,
        mut sequence: Vec<u8>,
        contexts: &SnpContexts,
    ) -> Result<DecodedContig, String> {
        let idx = color as usize;
        // Edited SNPs are put back on the stored strand.
        if let Some(snp_edits) = &self.snp_edits {
            snp_edits.revert(color, &mut sequence, &|f, pos| contexts.get(&(f, pos)).copied()).map_err(|e| e.to_string())?;
        }
        if let Some(runs) = self.n_runs.as_ref().map(|r| &r.runs[idx]).filter(|r| !r.is_empty()) {
            sequence = nruns::restore(&sequence, runs).ok_or_else(|| format!("the runs of N of colour {} do not fit it", color))?;
        }
        // Contigs stored on the other strand are reverse complemented back.
        let sequence = if self.reversed.as_ref().is_some_and(|r| r[idx]) {
            crate::strands::reverse_complement(&sequence)
        } else {
            sequence
        };
        Ok(DecodedContig { color, sequence })
    }

    /// The uncollapsed graph and the colours using each of its edges
    /// ([Csr::colorsets], recovered by walking every colour), for viewing.
    pub fn into_colored_graph(
        self,
    ) -> Result<(Csr<u32>, Vec<Vec<u32>>), E> {
        let csr = self.graph.to_csr();
        let colorsets = csr.colorsets(&self.choices, self.n_contigs)?;
        Ok((csr, colorsets))
    }
}

/// Decodes chosen contigs of an archive in chunks, see
/// [ArchiveReader::decoder]. Keeps the contexts that repeated SNPs of later
/// chunks need from the contigs decoded so far.
pub struct ContigDecoder<'r, 'a> {
    reader: &'r ArchiveReader<'a>,
    wanted: HashMap<u32, Vec<u32>>,
    contexts: SnpContexts,
    /// Per input file not complete yet: the digests of its contigs decoded
    /// so far.
    pending: HashMap<usize, Vec<Option<ContigDigest>>>,
    /// The checks of the files completed since [ContigDecoder::finished_files].
    finished: Vec<FileCheck>,
}

impl ContigDecoder<'_, '_> {
    /// Keeps the contexts later colours need from the edited sequence of
    /// `color`.
    fn record(
        &mut self,
        color: u32,
        edited: &[u8],
    ) {
        for &pos in self.wanted.get(&color).map_or(&[][..], |v| v.as_slice()) {
            if let Some(context) = snps::context_at(edited, pos as usize) {
                self.contexts.insert((color, pos), context);
            }
        }
    }

    /// Decodes the contigs `colors` in parallel, also decoding the contigs
    /// (not decoded yet) whose contexts their repeated SNPs need. An error
    /// if such a contig cannot be decoded; otherwise one result per colour.
    ///
    /// A file is checked against its hash once all its contigs have been
    /// decoded, over any calls and in any order ([ContigDecoder::finished_files]);
    /// the contigs of a file decoded in part are not checked.
    pub fn decode(
        &mut self,
        colors: &[u32],
    ) -> Result<Vec<Result<DecodedContig, ArchiveError>>, E> {
        let reader = self.reader;
        let edited: Vec<Result<Vec<u8>, String>> = colors.par_iter().map(|&color| reader.decode_edited(color)).collect();
        for (&color, sequence) in colors.iter().zip(&edited) {
            if let Ok(sequence) = sequence {
                self.record(color, sequence);
            }
        }
        // Contexts from colours not decoded yet (random access).
        if let Some(edits) = &reader.snp_edits {
            let mut missing: Vec<u32> = colors.iter().flat_map(|&c| edits.context_sources(c)).filter(|key| !self.contexts.contains_key(key)).map(|(f, _)| f).collect();
            missing.sort_unstable();
            missing.dedup();
            let sources: Vec<(u32, Result<Vec<u8>, String>)> = missing.par_iter().map(|&f| (f, reader.decode_edited(f))).collect();
            for (f, sequence) in sources {
                let sequence = sequence.map_err(|e| ArchiveError(format!("decoding color {}: {}", f, e)))?;
                self.record(f, &sequence);
            }
        }
        let contexts = &self.contexts;
        let checked = reader.hashes.is_some();
        let decoded: Vec<(Result<DecodedContig, ArchiveError>, Option<ContigDigest>)> = colors.par_iter().zip(edited.into_par_iter()).map(|(&color, sequence)| {
            let contig = sequence.and_then(|sequence| reader.finish_color(color, sequence, contexts))
                .map_err(|e| ArchiveError(format!("decoding color {}: {}", color, e)));
            let digest = contig.as_ref().ok().filter(|_| checked).map(|c| contig_digest(reader.contig_name(color).unwrap_or_default(), &c.sequence));
            (contig, digest)
        }).collect();
        for (&color, (_, digest)) in colors.iter().zip(&decoded) {
            if let Some(digest) = digest {
                self.add_digest(color, *digest);
            }
        }
        Ok(decoded.into_iter().map(|(contig, _)| contig).collect())
    }

    /// Records the digest of decoded contig `color`; checks its file if that
    /// completes it.
    fn add_digest(
        &mut self,
        color: u32,
        digest: ContigDigest,
    ) {
        let reader = self.reader;
        let Some(hashes) = &reader.hashes else { return };
        let file = reader.file_of(color);
        let range = reader.file_colors(file);
        let digests = self.pending.entry(file).or_insert_with(|| vec![None; range.len()]);
        digests[(color - range.start) as usize] = Some(digest);
        if digests.iter().all(|d| d.is_some()) {
            let digests: Vec<ContigDigest> = self.pending.remove(&file).unwrap().into_iter().flatten().collect();
            let hash_match = hashes.get(file) == Some(&file_hash(&digests));
            self.finished.push(FileCheck { file, hash_match });
        }
    }

    /// The checks of the files completed since the last call, in the order
    /// they were completed.
    pub fn finished_files(
        &mut self,
    ) -> Vec<FileCheck> {
        std::mem::take(&mut self.finished)
    }
}
