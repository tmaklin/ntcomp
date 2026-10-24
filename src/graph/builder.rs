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
//! Building the graph from the contigs' NLZ dictionaries.
//!
use core::ops::Range;

use std::collections::HashMap;

use indexmap::IndexSet;

use rayon::prelude::*;

use crate::edits::{self, CalledEdit};
use crate::strands::{self, StrandChooser};
use crate::index::Index;

use super::choices::{ChoiceWriter, Choices};
use super::collapsed::ChainCsr;
use super::snps::{self, SnpEdits};
use crate::nlz::{self, IMPLICIT_SUFFIX_LEN};
use super::{Csr, GraphWalk, StoredGraph, E, SINK_COLEX, SINK_NODE, SOURCE_COLEX, SOURCE_NODE};

/// Builds the graph and remembers the path each colour takes.
pub struct GraphBuilder {

    /// colex rank -> node index (position in the set)
    nodes: IndexSet<u32>,

    /// (source node, target node, suffix length); edge id is the position in the set
    edges: IndexSet<(u32, u32, u32)>,

    /// Per edge id: true if the suffix length is implicit, see [SuffixLen](crate::nlz::SuffixLen)
    implicit: Vec<bool>,

    /// Per colour: edge ids (indexes into `edges`) in traversal order
    paths: Vec<Vec<usize>>,

    /// Store linear chains as single edges in [GraphBuilder::finish]
    collapse_chains: bool,

    /// Choose the strand of every contig added with [GraphBuilder::add_contig]
    reorient_contigs: bool,

    /// One in this many _k_-mers is sampled for choosing the strands, see [strands]
    reorient_sampling: u64,

    /// Created when the first contig is added with reorienting on
    strand_chooser: Option<StrandChooser>,

    /// Per colour: true if the contig is stored reverse complemented
    reversed: Vec<bool>,

    /// Edit SNPs towards the index in [GraphBuilder::add_contig]
    edit_snps: bool,

    /// Edit SNPs shared by at most this many contigs, see [snps]
    snp_max_shared: usize,

    /// Per colour from [GraphBuilder::call_edits]: the chosen strand and the
    /// SNPs on it, until the contig is added
    called: Vec<(bool, Vec<CalledEdit>)>,

    /// Number of contigs with each SNP (by [CalledEdit::key])
    snp_counts: HashMap<u64, usize>,

    /// Per colour: the edited SNPs
    edited: Vec<Vec<snps::EditedSnp>>,
}

impl Default for GraphBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl GraphBuilder {
    pub fn new() -> Self {
        let mut nodes = IndexSet::new();
        nodes.insert(SOURCE_COLEX);
        nodes.insert(SINK_COLEX);

        GraphBuilder {
            nodes,
            edges: IndexSet::new(),
            implicit: Vec::new(),
            paths: Vec::new(),
            collapse_chains: false,
            reorient_contigs: false,
            reorient_sampling: strands::DEFAULT_SAMPLING,
            strand_chooser: None,
            reversed: Vec::new(),
            edit_snps: false,
            snp_max_shared: snps::DEFAULT_MAX_SHARED,
            called: Vec::new(),
            snp_counts: HashMap::new(),
            edited: Vec::new(),
        }
    }

    /// Collapse linear chains into single edges in [GraphBuilder::finish].
    pub fn collapse_chains(
        mut self,
        collapse: bool,
    ) -> Self {
        self.collapse_chains = collapse;
        self
    }

    /// Store every contig added with [GraphBuilder::add_contig] on the strand
    /// that shares the most with the contigs before it, see [strands].
    pub fn reorient_contigs(
        mut self,
        reorient: bool,
    ) -> Self {
        self.reorient_contigs = reorient;
        self
    }

    /// Sample one in `one_in` _k_-mers for choosing the strands
    /// ([strands::DEFAULT_SAMPLING] by default); panics if `one_in` is 0.
    pub fn reorient_sampling(
        mut self,
        one_in: u64,
    ) -> Self {
        assert!(one_in > 0, "sampling width must be at least 1");
        self.reorient_sampling = one_in;
        self
    }

    /// Edit SNPs towards the index before encoding, see [snps].
    ///
    /// Needs two passes over the contigs: [GraphBuilder::call_edits] for every
    /// contig first, then [GraphBuilder::add_contig] for every contig in the
    /// same order. The edited SNPs, [GraphBuilder::snp_edits], have to be put
    /// back after decoding.
    pub fn edit_snps(
        mut self,
        edit: bool,
    ) -> Self {
        self.edit_snps = edit;
        self
    }

    /// Edit only the SNPs that at most `max_shared` contigs have
    /// ([snps::DEFAULT_MAX_SHARED] by default); panics if `max_shared` is 0.
    pub fn edit_snps_max_shared(
        mut self,
        max_shared: usize,
    ) -> Self {
        assert!(max_shared > 0, "at least one contig must be able to have a SNP");
        self.snp_max_shared = max_shared;
        self
    }

    /// Per colour: true if the contig is stored reverse complemented.
    ///
    /// Only contigs added with [GraphBuilder::add_contig] with
    /// [GraphBuilder::reorient_contigs] on can be reversed. Their decoded
    /// sequences have to be reverse complemented to get the input back.
    pub fn reversed(&self) -> &[bool] {
        &self.reversed
    }

    /// Edited SNPs of every colour, None if no SNP was edited.
    ///
    /// They are put back into the decoded contig with [SnpEdits::revert]
    /// before reverse complementing it.
    pub fn snp_edits(&self) -> Option<SnpEdits> {
        if self.edited.iter().all(|snps| snps.is_empty()) {
            return None
        }
        Some(SnpEdits::with_repeats(&self.edited))
    }

    /// Chooses the strands of `sequences`, in order, if
    /// [GraphBuilder::reorient_contigs] is on. The _k_-mers are sampled in
    /// parallel.
    fn choose_strands(
        &mut self,
        sequences: &[&[u8]],
    ) -> Vec<bool> {
        if !self.reorient_contigs {
            return vec![false; sequences.len()]
        }
        let sampling = self.reorient_sampling;
        let chooser = self.strand_chooser.get_or_insert_with(|| StrandChooser::new(sampling));
        let sketches: Vec<_> = sequences.par_iter().map(|sequence| chooser.sketch(sequence)).collect();
        sketches.into_iter().map(|sketch| chooser.choose_sketched(sketch)).collect()
    }

    /// First pass of [GraphBuilder::edit_snps]: chooses the strand of the
    /// next contig and calls its SNPs on that strand. Returns the colour the
    /// contig will get from [GraphBuilder::add_contig].
    ///
    /// Arguments as in [GraphBuilder::add_contig].
    pub fn call_edits(
        &mut self,
        sequence: &[u8],
        index: &Index,
    ) -> u32 {
        self.call_edits_batch(&[sequence], index)[0]
    }

    /// Same as [GraphBuilder::call_edits] for the next `sequences.len()`
    /// contigs, in order. The SNPs are called in parallel; the result is the
    /// same as calling [GraphBuilder::call_edits] for each.
    pub fn call_edits_batch(
        &mut self,
        sequences: &[&[u8]],
        index: &Index,
    ) -> Vec<u32> {
        let strands = self.choose_strands(sequences);
        let called: Vec<Vec<CalledEdit>> = sequences.par_iter().zip(&strands).map(|(sequence, &reversed)| {
            if reversed {
                edits::call_edits(&strands::reverse_complement(sequence), index)
            } else {
                edits::call_edits(sequence, index)
            }
        }).collect();
        let first: u32 = self.called.len().try_into().unwrap();
        for (snps, reversed) in called.into_iter().zip(strands) {
            for snp in &snps {
                *self.snp_counts.entry(snp.key).or_default() += 1;
            }
            self.called.push((reversed, snps));
        }
        (first..(first + sequences.len() as u32)).collect()
    }

    /// Encodes the next contig and adds its path.
    ///
    /// With [GraphBuilder::reorient_contigs] on, chooses the strand of the contig
    /// first and encodes that strand. With [GraphBuilder::edit_snps] on, the
    /// strand comes from [GraphBuilder::call_edits], which must have been called
    /// for this contig, and the SNPs that at most
    /// [GraphBuilder::edit_snps_max_shared] contigs have are edited before
    /// encoding. `sequence` is a normalized nucleotide sequence, and `index`
    /// is the index passed to [nlz::parse].
    pub fn add_contig(
        &mut self,
        sequence: &[u8],
        index: &Index,
    ) -> Result<u32, E> {
        Ok(self.add_contigs(&[sequence], index)?[0])
    }

    /// Same as [GraphBuilder::add_contig] for the next `sequences.len()`
    /// contigs, in order. The contigs are encoded in parallel; the graph is
    /// the same as from adding them one at a time.
    pub fn add_contigs(
        &mut self,
        sequences: &[&[u8]],
        index: &Index,
    ) -> Result<Vec<u32>, E> {
        let first_color = self.paths.len();
        // Strands and SNPs to edit, in order.
        let prepared: Vec<(bool, Vec<CalledEdit>)> = if self.edit_snps {
            // The strands were chosen in the first pass.
            self.strand_chooser = None;
            if self.called.len() < first_color + sequences.len() {
                return Err(format!("SNPs of colour {} were not called before adding it", self.called.len().max(first_color)).into())
            }
            let (counts, max_shared) = (&self.snp_counts, self.snp_max_shared);
            self.called[first_color..(first_color + sequences.len())].iter_mut().map(|(reversed, called)| {
                let edited = std::mem::take(called).into_iter()
                    .filter(|snp| counts[&snp.key] <= max_shared)
                    .collect();
                (*reversed, edited)
            }).collect()
        } else {
            self.choose_strands(sequences).into_iter().map(|reversed| (reversed, Vec::new())).collect()
        };

        let strands: Vec<bool> = prepared.iter().map(|(reversed, _)| *reversed).collect();
        type Encoded = (usize, Vec<(usize, Range<usize>)>, Vec<bool>, Vec<snps::EditedSnp>);
        let encoded: Vec<Result<Encoded, String>> = sequences.par_iter().zip(prepared.into_par_iter()).enumerate().map(|(i, (sequence, (reversed, edited)))| {
            let oriented = if reversed { strands::reverse_complement(sequence) } else { sequence.to_vec() };
            let edits = edited.into_iter().map(|called| called.edit).collect();
            let (oriented, edits) = edits::apply_edits(&oriented, edits)
                .ok_or_else(|| format!("contig {} differs from the one its SNPs were called on", first_color + i))?;
            // Contexts in the edited contig, for storing shared SNPs once.
            let edited = snps::with_contexts(&oriented, edits);
            let (dictionary, implicit) = nlz::parse(&oriented, index).map_err(|e| format!("contig {}: {}", first_color + i, e))?;
            Ok((oriented.len(), dictionary, implicit, edited))
        }).collect();

        let mut colors = Vec::with_capacity(sequences.len());
        for (encoded, reversed) in encoded.into_iter().zip(strands) {
            let (len, dictionary, implicit, edited) = encoded?;
            let color = self.add_path_with_flags(len, &dictionary, &implicit);
            self.reversed[color as usize] = reversed;
            self.edited[color as usize] = edited;
            colors.push(color);
        }
        Ok(colors)
    }

    fn insert_node(
        &mut self,
        colex_rank: usize
    ) -> u32 {
        let colex_rank: u32 = colex_rank.try_into().expect("colex rank does not fit in u32");
        assert!(colex_rank < SOURCE_COLEX, "colex rank collides with dummy node ranks");
        self.nodes.insert_full(colex_rank).0 as u32
    }

    fn edge(
        &mut self,
        source: u32,
        target: u32,
        weight: u32,
        implicit: bool,
    ) -> usize {
        let (id, inserted) = self.edges.insert_full((source, target, weight));
        if inserted {
            self.implicit.push(false);
        }
        // Implicit if any colour saw the extension stop in the SBWT: the length
        // then only depends on the _k_-mer, not on where the colour uses it.
        self.implicit[id] |= implicit;
        id
    }

    /// Adds the dictionary of the next colour.
    ///
    /// `dictionary` is the output of [nlz::parse]. The edge out of
    /// the source node stores `sequence_length`.
    pub fn add_path(
        &mut self,
        sequence_length: usize,
        dictionary: &[(usize, Range<usize>)],
    ) -> u32 {
        self.add_path_with_flags(sequence_length, dictionary, &vec![false; dictionary.len()])
    }

    /// Same as [GraphBuilder::add_path] but stores the suffix length of the
    /// entries flagged in `implicit` as [SuffixLen::Implicit](crate::nlz::SuffixLen::Implicit).
    ///
    /// `dictionary` and `implicit` are the output of [nlz::parse].
    pub fn add_path_with_flags(
        &mut self,
        sequence_length: usize,
        dictionary: &[(usize, Range<usize>)],
        implicit: &[bool],
    ) -> u32 {
        assert_eq!(dictionary.len(), implicit.len());
        let color: u32 = self.paths.len().try_into().unwrap();
        let mut path = Vec::with_capacity(dictionary.len() + 1);

        let mut prev = SOURCE_NODE;
        let mut prev_weight: u32 = sequence_length.try_into().expect("sequence length does not fit in u32");
        let mut prev_implicit = false;
        for ((suffix_len, colex), &is_implicit) in dictionary.iter().zip(implicit) {
            let node = self.insert_node(colex.start);
            path.push(self.edge(prev, node, prev_weight, prev_implicit));
            prev = node;
            prev_weight = (*suffix_len).try_into().unwrap();
            assert!(prev_weight != 0, "empty dictionary entry");
            prev_implicit = is_implicit;
        }
        path.push(self.edge(prev, SINK_NODE, prev_weight, prev_implicit));

        self.paths.push(path);
        self.reversed.push(false);
        self.edited.push(Vec::new());
        color
    }

    /// Converts into the stored representation (graph + choice stream),
    /// collapsing linear chains if [GraphBuilder::collapse_chains] was set.
    pub fn finish(self) -> (StoredGraph, Choices) {
        if self.collapse_chains {
            let (chains, choices) = self.to_chains();
            (StoredGraph::Chains(chains), choices)
        } else {
            let (csr, choices) = self.to_csr();
            (StoredGraph::Csr(csr), choices)
        }
    }

    /// Converts into the stored representation (graph + choice stream).
    pub fn to_csr(self) -> (Csr<u32>, Choices) {
        let (csr, csr_pos) = self.build_csr();
        let choices = self.choices(&csr, &csr_pos);
        (csr, choices)
    }

    /// Converts into the stored representation with linear chains collapsed.
    ///
    /// Only junctions have more than one outgoing edge and their chains keep
    /// the CSR order of their first edges, so the choice stream is the same as
    /// for [GraphBuilder::to_csr].
    pub fn to_chains(self) -> (ChainCsr, Choices) {
        let (csr, csr_pos) = self.build_csr();
        let choices = self.choices(&csr, &csr_pos);
        (ChainCsr::from_csr(&csr), choices)
    }

    /// Records which edge each colour takes at every node with more than one outgoing edge.
    fn choices(
        &self,
        csr: &Csr<u32>,
        csr_pos: &[u32],
    ) -> Choices {
        let streams = self.paths.iter().map(|path| {
            let mut writer = ChoiceWriter::new();
            for &id in path {
                let source = self.edges.get_index(id).unwrap().0 as usize;
                let (start, end) = csr.row(source);
                if end - start > 1 {
                    writer.push((csr_pos[id] as usize - start) as u32, (end - start) as u32);
                }
            }
            writer.finish()
        }).collect();
        Choices { streams }
    }

    /// Builds the CSR and the CSR position of every edge id.
    fn build_csr(&self) -> (Csr<u32>, Vec<u32>) {
        let n_nodes = self.nodes.len();
        let n_edges = self.edges.len();

        // Number of times the colours take each edge.
        let mut uses = vec![0_u32; n_edges];
        for path in &self.paths {
            for &id in path {
                uses[id] += 1;
            }
        }

        // Edge order: by source, then the most used edges first, then target
        // and suffix length. The choice streams store ranks in this order, and
        // rank 0 being the most common choice makes them compress well.
        let mut order: Vec<usize> = (0..n_edges).collect();
        order.sort_unstable_by_key(|&id| {
            let &(source, target, weight) = self.edges.get_index(id).unwrap();
            (source, std::cmp::Reverse(uses[id]), target, weight)
        });

        let mut csr_pos = vec![0_u32; n_edges];
        let mut row_ptrs = vec![0_u32; n_nodes + 1];
        let mut column_indices = Vec::with_capacity(n_edges);
        let mut suffix_lens = Vec::with_capacity(n_edges);

        for (pos, &id) in order.iter().enumerate() {
            let &(source, target, weight) = self.edges.get_index(id).unwrap();
            csr_pos[id] = pos as u32;
            row_ptrs[source as usize + 1] += 1;
            column_indices.push(target);
            suffix_lens.push(if self.implicit[id] { IMPLICIT_SUFFIX_LEN } else { weight });
        }
        for i in 0..n_nodes {
            row_ptrs[i + 1] += row_ptrs[i];
        }

        let csr = Csr {
            node_weights: self.nodes.iter().copied().collect(),
            row_ptrs,
            column_indices,
            suffix_lens,
        };

        (csr, csr_pos)
    }
}

