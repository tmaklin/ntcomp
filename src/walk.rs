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
//! Walks through the de Bruijn graph of the index. (Walks along a colour's
//! path through the stored graph are [GraphWalk](crate::graph::GraphWalk).)
//!
//! [IndexWalk] walks the de Bruijn graph of an SBWT index. Its operations are
//! default methods built on two accessors, the SBWT and the LCS value of a
//! row, so they are shared by every way of reaching the index:
//! [Index](crate::index::Index) (the SBWT with its LCS array) and [WalkIndex] (the LCS array
//! optional, the values otherwise computed from the SBWT with
//! [lcs_from_sbwt]). Used by the parser ([crate::nlz::parse]), the
//! decoder ([crate::nlz::decode_dictionary]), the node walks of the stored graph
//! ([crate::graph::collapsed]) and the edit caller
//! ([crate::edits::call_edits]).
//!
//! The k-mers of a _group_ share their last _k_ − 1 characters: they are
//! adjacent in colex order, with LCS values of at least _k_ − 1 between them,
//! and are the predecessors of the same _k_-mers. The inverse LF step of a
//! row gives the first row of its predecessors' group.

use std::ops::Range;

use sbwt::{ExtendRight, LcsArray, SbwtIndex, SubsetMatrix, SubsetSeq};

/// A pattern in the index while stepping right: its colex interval and its
/// length (at most _k_).
pub type Pattern = (Range<usize>, usize);

/// Walks through the de Bruijn graph of an SBWT index (see the module docs).
/// The operations need the select support of the SBWT.
pub trait IndexWalk {
    /// The SBWT.
    fn sbwt(&self) -> &SbwtIndex<SubsetMatrix>;

    /// The LCS array value of `row`: the longest common suffix of the
    /// _k_-mers at `row` − 1 and `row`.
    fn lcs_at(&self, row: usize) -> usize;

    /// The _k_-mer size.
    fn k(&self) -> usize {
        self.sbwt().k()
    }

    /// `interval` widened to whole groups: the rows sharing the last _k_ − 1
    /// characters with its first and last row.
    fn widen(
        &self,
        interval: Range<usize>,
    ) -> Range<usize> {
        let k = self.k();
        let (start, mut end) = (self.group_start(interval.start), interval.end);
        while end < self.sbwt().n_sets() && self.lcs_at(end) >= k - 1 {
            end += 1;
        }
        start..end
    }

    /// The first row of the group of `row`: the one that stores the group's
    /// outgoing edges.
    fn group_start(
        &self,
        row: usize,
    ) -> usize {
        let k = self.k();
        let mut start = row;
        while start > 0 && self.lcs_at(start) >= k - 1 {
            start -= 1;
        }
        start
    }

    /// The group of `row`: the rows sharing its last _k_ − 1 characters.
    fn group(
        &self,
        row: usize,
    ) -> Range<usize> {
        self.widen(row..row + 1)
    }

    /// Whether `row`, the first of its group (as from an inverse LF step), is
    /// alone in it: the only predecessor of the _k_-mers it precedes.
    fn only_in_group(
        &self,
        row: usize,
    ) -> bool {
        let k = self.k();
        (row == 0 || self.lcs_at(row) < k - 1) && (row + 1 >= self.sbwt().n_sets() || self.lcs_at(row + 1) < k - 1)
    }

    /// The predecessors of the _k_-mer at `row` (one group); None for the
    /// source row.
    fn predecessors(
        &self,
        row: usize,
    ) -> Option<Range<usize>> {
        Some(self.group(self.sbwt().inverse_lf_step(row)?))
    }

    /// Number of successors of the _k_-mer at `row` (the edges stored at its
    /// group's first row).
    fn out_degree(
        &self,
        row: usize,
    ) -> usize {
        self.sbwt().sbwt().subset_size(self.group_start(row))
    }

    /// Pushes every successor of the _k_-mer at `row` with the base that leads
    /// to it, in A, C, G, T order, to `out`.
    fn push_successors(
        &self,
        row: usize,
        out: &mut Vec<(usize, u8)>,
    ) {
        let sbwt = self.sbwt();
        let start = self.group_start(row);
        for (i, &c) in sbwt.alphabet().iter().enumerate() {
            if sbwt.sbwt().set_contains(start, i as u8) {
                out.push((sbwt.lf_step(start, i), c));
            }
        }
    }

    /// Number of real predecessors of the _k_-mer at `row` (dummy rows padded
    /// with $ not counted). A _k_-mer with no real predecessor has one dummy
    /// predecessor, alone in its group, so only a single-row group is checked.
    fn in_degree(
        &self,
        row: usize,
    ) -> usize {
        match self.predecessors(row) {
            None => 0,
            Some(group) if group.len() == 1 => self.is_real_kmer(group.start) as usize,
            Some(group) => group.len(),
        }
    }

    /// Whether the _k_-mer at `row` is a real _k_-mer, not a dummy one padded
    /// with $: the _k_ − 1 inverse LF steps from it do not reach the source
    /// row.
    fn is_real_kmer(
        &self,
        row: usize,
    ) -> bool {
        let sbwt = self.sbwt();
        let mut cur = row;
        for _ in 0..sbwt.k() {
            match sbwt.inlabel(cur).and_then(|_| sbwt.inverse_lf_step(cur)) {
                Some(predecessor) => cur = predecessor,
                None => return false,
            }
        }
        true
    }

    /// Left extension of the _k_-mer at `row`: the labels read along inverse
    /// LF steps from `row` are the _k_-mer's bases from right to left,
    /// followed by the first bases of the predecessors taken, so a phrase of
    /// length _k_ + L takes _k_ + L steps and no searches.
    ///
    /// With `check`, extends while the _k_-mer has exactly one predecessor, at
    /// most `max_extension_len` characters, and returns whether the extension
    /// stopped there: the same as [crate::nlz::left_extend_kmer] without a
    /// reference on the _k_-mer at `row`. Without `check`, takes
    /// `max_extension_len` steps without looking at the in-degrees, for a
    /// length known to be an extension (each _k_-mer on the way has one
    /// predecessor); the returned flag is then false.
    ///
    /// Panics without `check` if the walk reaches a dummy _k_-mer.
    fn left_extend(
        &self,
        row: usize,
        max_extension_len: usize,
        check: bool,
    ) -> (Vec<u8>, bool) {
        let sbwt = self.sbwt();
        let k = sbwt.k();
        // labels[t] is the label of the row t steps from `row`: the base t from the
        // right of the phrase. The row taken at extension step j has its first base
        // at labels[j + k - 1].
        let mut labels: Vec<Option<u8>> = Vec::with_capacity(k + 64);
        labels.push(sbwt.inlabel(row));
        let mut cur = Some(row);
        let mut extension_len = 0;
        let mut stopped_in_sbwt = false;
        while extension_len < max_extension_len {
            let Some(predecessor) = cur.and_then(|r| sbwt.inverse_lf_step(r)) else {
                stopped_in_sbwt = true;
                break
            };
            if check && !self.only_in_group(predecessor) {
                stopped_in_sbwt = true;
                break
            }
            labels.push(sbwt.inlabel(predecessor));
            cur = Some(predecessor);
            extension_len += 1;
        }
        while labels.len() < k + extension_len {
            cur = cur.and_then(|r| sbwt.inverse_lf_step(r));
            labels.push(cur.and_then(|r| sbwt.inlabel(r)));
        }
        // A dummy predecessor (the k-mer has none in the index) has a $ as its
        // first base: the walk reaches the source row within k - 1 steps after it.
        if let Some(dollar) = labels[..(k + extension_len)].iter().position(|label| label.is_none()) {
            assert!(check && dollar >= k, "the walk from colex rank {} reached a dummy k-mer", row);
            extension_len = dollar - k;
            stopped_in_sbwt = true;
        }
        // Every label before the first $ is there.
        let phrase = labels[..(k + extension_len)].iter().rev().flatten().copied().collect();
        (phrase, stopped_in_sbwt)
    }

    /// The first _k_ bases of the phrase of the _k_-mer at `colex` with suffix
    /// length `len` (0: implicit, extended as far as it goes); None for a
    /// short entry (`len` < _k_) or a rank outside the index.
    fn phrase_start(
        &self,
        colex: u32,
        len: u32,
    ) -> Option<Vec<u8>> {
        let k = self.k();
        if (len != 0 && (len as usize) < k) || colex as usize >= self.sbwt().n_sets() {
            return None
        }
        let max = if len == 0 { usize::MAX / 2 } else { len as usize - k };
        let (phrase, _) = self.left_extend(colex as usize, max, len == 0);
        Some(phrase[..k].to_vec())
    }

    /// One step right: the interval of `pattern` followed by `c`, dropping the
    /// pattern's first character first if it is _k_ long. None if that is not
    /// in the index.
    fn step_right(
        &self,
        pattern: &Pattern,
        c: u8,
    ) -> Option<Pattern> {
        let k = self.k();
        let (mut interval, mut len) = pattern.clone();
        if len == k {
            interval = self.widen(interval);
            len = k - 1;
        }
        let next = self.sbwt().extend_right(interval, c);
        (!next.is_empty()).then_some((next, len + 1))
    }

    /// Steps right by the characters of `sequence` from `from`, stopping
    /// before `until` or at the first one not in the index. Returns the
    /// pattern and the first position not taken.
    fn walk_right(
        &self,
        sequence: &[u8],
        mut pattern: Pattern,
        from: usize,
        until: usize,
    ) -> (Pattern, usize) {
        let mut pos = from;
        while pos < until {
            match self.step_right(&pattern, sequence[pos]) {
                Some(next) => {
                    pattern = next;
                    pos += 1;
                },
                None => break,
            }
        }
        (pattern, pos)
    }
}

/// The SBWT with its LCS array if there is one; otherwise the LCS values are
/// computed from the SBWT ([lcs_from_sbwt], one step per common character).
/// Storing a graph with derived nodes needs the array (for the matching
/// statistics of the walks); reading one does not.
#[derive(Clone, Copy)]
pub struct WalkIndex<'a> {
    pub sbwt: &'a SbwtIndex<SubsetMatrix>,
    pub lcs: Option<&'a LcsArray>,
}

impl IndexWalk for WalkIndex<'_> {
    fn sbwt(&self) -> &SbwtIndex<SubsetMatrix> {
        self.sbwt
    }

    fn lcs_at(&self, row: usize) -> usize {
        match self.lcs {
            Some(lcs) => lcs.access(row),
            None => lcs_from_sbwt(self.sbwt, row),
        }
    }
}

/// LCS array value of `colex_rank` computed from the SBWT: the length of the
/// longest common suffix of the _k_-mers at `colex_rank` - 1 and
/// `colex_rank` (0 for rank 0). Takes one step per common character.
///
/// Requires the select support of the SBWT.
pub fn lcs_from_sbwt(
    sbwt: &SbwtIndex<SubsetMatrix>,
    colex_rank: usize,
) -> usize {
    if colex_rank == 0 {
        return 0
    }
    let (mut a, mut b) = (colex_rank - 1, colex_rank);
    for len in 0..sbwt.k() {
        if sbwt.inlabel(a).is_none() || sbwt.inlabel(a) != sbwt.inlabel(b) {
            return len
        }
        // Both rows have the same label, so neither is the source row.
        match (sbwt.inverse_lf_step(a), sbwt.inverse_lf_step(b)) {
            (Some(next_a), Some(next_b)) => (a, b) = (next_a, next_b),
            _ => return len,
        }
    }
    sbwt.k()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::related_index;

    #[test]
    fn lcs_array_and_lcs_from_sbwt_agree() {
        for k in [15, 31] {
            let (_, index) = related_index(k);
            let (sbwt, lcs) = (index.sbwt(), index.lcs());
            assert_eq!(lcs.len(), sbwt.n_sets());
            for rank in 0..lcs.len() {
                assert_eq!(lcs_from_sbwt(sbwt, rank), lcs.access(rank), "k {} rank {}", k, rank);
            }
        }
    }
}
