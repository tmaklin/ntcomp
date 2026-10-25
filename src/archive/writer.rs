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

use crate::edits;
use crate::fastx::Contig;
use crate::graph::encoder::GraphEncoder;
use crate::graph::snps::{self, EditedSnp, SnpEdits};
use crate::graph::{contig_hash, ContigHash, GraphBuilder};
use crate::index::Index;
use crate::nlz;
use crate::nruns::{self, NRun, NRuns};
use crate::strands::{self, StrandChooser};

type E = Box<dyn std::error::Error>;

/// An input file: its name and its contigs.
pub type InputFile = (Vec<u8>, Vec<Contig>);

/// Options of an archive. The defaults are the best archives: reorienting
/// and editing on, as in `ntcomp graph`.
#[derive(Clone, Debug, PartialEq)]
pub struct ArchiveOpts {
    /// Store every contig on the strand that shares the most with the
    /// contigs before it, see [strands].
    pub reorient: bool,
    /// One in this many _k_-mers is sampled for choosing the strands (at
    /// least 1).
    pub reorient_sampling: u64,
    /// Edit SNPs and small indels towards the index before encoding, see
    /// [edits].
    pub edit: bool,
}

impl Default for ArchiveOpts {
    fn default() -> Self {
        ArchiveOpts { reorient: true, reorient_sampling: strands::DEFAULT_SAMPLING, edit: true }
    }
}

/// Builds an archive from input files added in order. Each contig gets the
/// next colour; the archive keeps its hash, its name and its file.
///
/// A contig is stored on the strand chosen for it ([ArchiveOpts::reorient]),
/// without its runs of `N` ([nruns]), with its edits applied
/// ([ArchiveOpts::edit]), as the path of its NLZ dictionary ([nlz::parse])
/// in the graph ([GraphBuilder]); the strands, the runs and the edits are
/// stored for putting the contig back after decoding. The pieces between
/// the runs are edited and parsed as separate sequences, so that nothing
/// spans a cut.
pub struct ArchiveWriter<'a> {
    opts: ArchiveOpts,
    index: &'a Index,
    builder: GraphBuilder,
    /// Created when the first strand is chosen.
    strand_chooser: Option<StrandChooser>,
    /// Per colour: whether the contig is stored reverse complemented.
    reversed: Vec<bool>,
    /// Per colour: the edits applied, with their contexts.
    edited: Vec<Vec<EditedSnp>>,
    /// Per colour: the runs of `N` cut out.
    n_runs: Vec<Vec<NRun>>,
    hashes: Vec<ContigHash>,
    contig_names: Vec<Vec<u8>>,
    file_ranges: Vec<(Vec<u8>, Range<u32>)>,
}

impl<'a> ArchiveWriter<'a> {
    /// A writer with the options `opts` that encodes the contigs against
    /// `index`.
    pub fn new(
        opts: ArchiveOpts,
        index: &'a Index,
    ) -> Self {
        assert!(opts.reorient_sampling > 0, "sampling width must be at least 1");
        let builder = GraphBuilder::new();
        ArchiveWriter {
            opts, index, builder,
            strand_chooser: None,
            reversed: Vec::new(),
            edited: Vec::new(),
            n_runs: Vec::new(),
            hashes: Vec::new(),
            contig_names: Vec::new(),
            file_ranges: Vec::new(),
        }
    }

    /// Chooses the strands of `sequences`, in order, if
    /// [ArchiveOpts::reorient] is on. The _k_-mers are sampled in parallel.
    fn choose_strands(
        &mut self,
        sequences: &[&[u8]],
    ) -> Vec<bool> {
        if !self.opts.reorient {
            return vec![false; sequences.len()]
        }
        let sampling = self.opts.reorient_sampling;
        let chooser = self.strand_chooser.get_or_insert_with(|| StrandChooser::new(sampling));
        let sketches: Vec<_> = sequences.par_iter().map(|sequence| chooser.sketch(sequence)).collect();
        sketches.into_iter().map(|sketch| chooser.choose_sketched(sketch)).collect()
    }

    /// Adds the contigs of `files`, encoded in parallel. Returns the colours
    /// they got.
    ///
    /// A contig's matching statistics are computed once: the edits are
    /// called from them ([ArchiveOpts::edit]), and after the edits they are
    /// updated from each edit on ([edits::matching_statistics_after_edits])
    /// for parsing the edited contig.
    pub fn add_files(
        &mut self,
        files: &[InputFile],
    ) -> Result<Vec<u32>, E> {
        let sequences: Vec<&[u8]> = files.iter().flat_map(|(_, contigs)| contigs).map(|(_, sequence)| sequence.as_slice()).collect();
        let first_color = self.reversed.len();
        let strands = self.choose_strands(&sequences);

        let (index, edit) = (self.index, self.opts.edit);
        type Encoded = (bool, nlz::Dictionary, Vec<EditedSnp>, Vec<NRun>);
        let encoded: Vec<Result<Encoded, String>> = sequences.par_iter().zip(strands.into_par_iter()).enumerate().map(|(i, (sequence, reversed))| {
            let color = first_color + i;
            let oriented = if reversed { strands::reverse_complement(sequence) } else { sequence.to_vec() };
            let (cut, runs) = nruns::cut(&oriented).ok_or_else(|| format!("contig {}: a character other than A, C, G, T, a gap or an IUPAC code, or a run of them too long", color))?;
            // Each piece between the runs is edited and parsed on its own;
            // the edited pieces, their edits and their dictionaries are
            // joined (the dictionary from the end: the last piece first).
            let mut edited: Vec<u8> = Vec::with_capacity(cut.len());
            let mut all_edits: Vec<edits::Edit> = Vec::new();
            let mut pieces: Vec<nlz::Dictionary> = Vec::new();
            for segment in nruns::segments(cut.len(), &runs) {
                let piece = &cut[segment];
                let ms = index.matching_statistics(piece);
                let edits = if edit { edits::call_edits_with(piece, &ms, index) } else { Vec::new() };
                let (edited_piece, edits) = edits::apply_edits(piece, edits)
                    .ok_or_else(|| format!("contig {}: the edits do not match it", color))?;
                let ms = edits::matching_statistics_after_edits(ms, &edited_piece, &edits, index);
                pieces.push(nlz::parse_with(&edited_piece, ms, index).map_err(|e| format!("contig {}: {}", color, e))?);
                let offset = u32::try_from(edited.len()).map_err(|_| format!("contig {}: longer than 2^32 bases", color))?;
                all_edits.extend(edits.into_iter().map(|e| edits::Edit { pos: e.pos + offset, ..e }));
                edited.extend_from_slice(&edited_piece);
            }
            let dictionary = nlz::Dictionary {
                sequence_len: edited.len(),
                entries: pieces.into_iter().rev().flat_map(|d| d.entries).collect(),
            };
            // Contexts in the edited contig, for storing shared SNPs once.
            let edited_snps = snps::with_contexts(&edited, all_edits);
            Ok((reversed, dictionary, edited_snps, runs))
        }).collect();

        let encoded = encoded.into_iter().collect::<Result<Vec<Encoded>, String>>()?;
        self.hashes.par_extend(sequences.par_iter().map(|sequence| contig_hash(sequence)));
        let mut colors = Vec::with_capacity(sequences.len());
        for (reversed, dictionary, edited, runs) in encoded {
            colors.push(self.builder.add_dictionary(&dictionary));
            self.reversed.push(reversed);
            self.edited.push(edited);
            self.n_runs.push(runs);
        }

        let mut color = self.contig_names.len() as u32;
        for (file_name, contigs) in files {
            let start = color;
            self.contig_names.extend(contigs.iter().map(|(name, _)| name.clone()));
            color += contigs.len() as u32;
            self.file_ranges.push((file_name.clone(), start..color));
        }
        Ok(colors)
    }

    /// Per colour added so far: whether the contig is stored reverse
    /// complemented.
    pub fn reversed(&self) -> &[bool] {
        &self.reversed
    }

    /// The edits of every colour as stored, None if no contig was edited.
    fn snp_edits(&self) -> Option<SnpEdits> {
        if self.edited.iter().all(|snps| snps.is_empty()) {
            return None
        }
        Some(SnpEdits::with_repeats(&self.edited))
    }

    /// The runs of `N` of every colour as stored, None if no contig has any.
    fn n_runs(&self) -> Option<NRuns> {
        if self.n_runs.iter().all(|runs| runs.is_empty()) {
            return None
        }
        Some(NRuns { runs: self.n_runs.clone() })
    }

    /// Writes the archive to `out`. Returns a summary of the section sizes.
    pub fn write<W: Write>(
        self,
        out: &mut W,
    ) -> Result<String, E> {
        let snp_edits = self.snp_edits();
        let n_runs = self.n_runs();
        let (graph, choices) = self.builder.finish();
        let mut encoder = GraphEncoder::new(graph, choices)
            .with_hashes(self.hashes)
            .with_contig_names(self.contig_names)
            .with_file_ranges(self.file_ranges)
            .with_reversed(self.reversed)
            .with_index(self.index);
        if let Some(snp_edits) = snp_edits {
            encoder = encoder.with_snp_edits(snp_edits);
        }
        if let Some(n_runs) = n_runs {
            encoder = encoder.with_n_runs(n_runs);
        }
        encoder.encode_to_write(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archive::ArchiveReader;
    use crate::edits::base_code;
    use crate::index::BuildOpts;
    use crate::testing::random_sequence;

    #[test]
    fn shared_edits_are_stored_once_and_put_back() {
        let reference = random_sequence(4000, 7);
        // SNPs are only called after matches longer than random ones.
        let index = Index::build(std::slice::from_ref(&reference), &BuildOpts { k: 63, ..Default::default() });

        // Two contigs with the same SNPs and one with its own, more than k apart.
        let with_snps = |positions: &[usize]| {
            let mut contig = reference.clone();
            for &i in positions {
                contig[i] = b"ACGT"[(base_code(contig[i]).unwrap() as usize + 1) % 4];
            }
            contig
        };
        let contigs: Vec<Contig> = [with_snps(&[500, 1500, 2500]), with_snps(&[500, 1500, 2500]), with_snps(&[1000, 3000])]
            .into_iter().enumerate().map(|(i, contig)| (format!("contig{}", i).into_bytes(), contig)).collect();
        let files: Vec<InputFile> = vec![(b"input".to_vec(), contigs.clone())];

        let opts = ArchiveOpts { reorient: false, ..Default::default() };
        let mut writer = ArchiveWriter::new(opts, &index);
        assert_eq!(writer.add_files(&files).unwrap(), vec![0, 1, 2]);
        let edits = writer.snp_edits().unwrap();
        assert_eq!(edits.n_snps(), 8);
        assert!(edits.n_repeats() > 0, "the second contig repeats the first one's SNPs");
        let mut bytes = Vec::new();
        writer.write(&mut bytes).unwrap();

        // In reverse order: the repeats need contexts from contigs not decoded yet.
        let reader = ArchiveReader::open(&mut bytes.as_slice(), &index).unwrap();
        let colors = [2, 1, 0];
        for decoded in reader.decoder(&colors).unwrap().decode(&colors).unwrap() {
            let decoded = decoded.unwrap();
            assert_eq!(decoded.sequence, contigs[decoded.color as usize].1, "colour {}", decoded.color);
            assert_eq!(decoded.hash_match, Some(true));
        }
    }
}
