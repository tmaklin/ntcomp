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
//! Editing contigs towards the index: SNPs and small indels.
//!
//! A contig that differs from the index by a SNP is encoded as several short
//! dictionary entries around the SNP. Replacing the contig's base with the
//! index's base before encoding gives fewer, longer entries. Small indels are
//! edited the same way: contig bases the index does not have are removed (a
//! deletion, [Change::Del]) and index bases the contig does not have are
//! added (an insertion, [Change::Ins]).
//!
//! The edits are called by walking the index ([call_edits]) and applied with
//! [apply_edits]. The archive stores them and puts them back after decoding
//! ([crate::graph::snps]); the tig builder leaves out gaps that editing
//! repairs ([crate::tigs]).
//!
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};


use crate::index::Index;
use crate::walk::{IndexWalk, Pattern};

/// Bases on both sides of an edit that identify it across contigs.
pub const FLANK: usize = 32;

/// A change of a contig against the index.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Change {
    /// A SNP: the contig's base `base` where the index has `index_base`.
    Sub { base: u8, index_base: u8 },
    /// Contig bases the index does not have, removed when editing.
    Del { bases: Vec<u8> },
    /// Index bases the contig does not have, added when editing.
    Ins { bases: Vec<u8> },
}

/// A change at a position of a contig: the first base it changes or removes,
/// or the base an insertion goes before.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Edit {
    pub pos: u32,
    pub change: Change,
}

/// A called edit (in the contig's coordinates) and the hash that identifies
/// it across contigs ([FLANK] contig bases on both sides and the change).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CalledEdit {
    pub edit: Edit,
    pub key: u64,
}

pub(crate) fn base_code(
    base: u8,
) -> Option<u8> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Hash of the edit with [FLANK] bases of `sequence` on both sides (before
/// its position and after the bases it changes or removes).
fn key(
    sequence: &[u8],
    edit: &Edit,
) -> u64 {
    let pos = edit.pos as usize;
    let after = match &edit.change {
        Change::Sub { .. } => pos + 1,
        Change::Del { bases } => pos + bases.len(),
        Change::Ins { .. } => pos,
    };
    let mut hasher = DefaultHasher::new();
    sequence[pos.saturating_sub(FLANK)..pos].hash(&mut hasher);
    edit.change.hash(&mut hasher);
    sequence[after..(after + FLANK).min(sequence.len())].hash(&mut hasher);
    hasher.finish()
}

/// Largest probability of a random match that edits are called after, see
/// [call_edits] (the default of kbo's variant caller).
const MAX_ERROR_PROB: f64 = 0.0000001;

/// Most edits in one stretch of the contig not in the index, see
/// [call_edits].
const MAX_EDITS: usize = 5;

/// Most contig bases one deletion removes, see [call_edits].
const MAX_DELETED: usize = 3;

/// Most index bases one insertion adds, see [call_edits].
const MAX_INSERTED: usize = 2;

/// An edit tried in a stretch the contig cannot extend: its change (as
/// the base, the number of contig bases removed, or the bases added), the
/// pattern after the walk following it, the contig position the walk
/// reached and whether that completes the stretch.
type Tried = (TriedChange, Pattern, usize, bool);

#[derive(Clone)]
enum TriedChange {
    Sub(u8),
    Del(usize),
    Ins(Vec<u8>),
}

/// Edits from `pos` until a walk takes _k_ contig bases from the last one,
/// starting from the match `pattern`: at most [MAX_EDITS], substitutions only
/// or also indels. None if the walk does not get through.
fn repair(
    index: &Index,
    sequence: &[u8],
    mut pattern: Pattern,
    mut pos: usize,
    with_indels: bool,
) -> Option<(Vec<(usize, TriedChange)>, Pattern, usize)> {
    let (k, n) = (index.k(), sequence.len());
    let mut found: Vec<(usize, TriedChange)> = Vec::new();
    loop {
        if found.len() == MAX_EDITS || pos + k > n || base_code(sequence[pos]).is_none() {
            return None
        }
        let mut tried: Vec<Tried> = Vec::new();
        for &b in b"ACGT" {
            if b == sequence[pos] {
                continue
            }
            let Some(after) = index.step_right(&pattern, b) else { continue };
            let (walked, reached) = index.walk_right(sequence, after, pos + 1, pos + k);
            tried.push((TriedChange::Sub(b), walked, reached, reached == pos + k));
        }
        if with_indels && !tried.iter().any(|t| t.3) {
            // Contig bases removed: the walk goes on after them.
            for d in 1..=MAX_DELETED {
                if pos + d + k > n || !sequence[pos..pos + d].iter().all(|&c| base_code(c).is_some()) {
                    break
                }
                let (walked, reached) = index.walk_right(sequence, pattern.clone(), pos + d, pos + d + k);
                if reached > pos + d {
                    tried.push((TriedChange::Del(d), walked, reached, reached == pos + d + k));
                }
            }
            // Index bases added: every path of up to MAX_INSERTED steps.
            let mut paths: Vec<(Vec<u8>, Pattern)> = vec![(Vec::new(), pattern.clone())];
            for _ in 0..MAX_INSERTED {
                let mut longer = Vec::new();
                for (bases, path) in &paths {
                    for &b in b"ACGT" {
                        let Some(after) = index.step_right(path, b) else { continue };
                        let mut bases = bases.clone();
                        bases.push(b);
                        let (walked, reached) = index.walk_right(sequence, after.clone(), pos, pos + k);
                        if reached > pos {
                            tried.push((TriedChange::Ins(bases.clone()), walked, reached, reached == pos + k));
                        }
                        longer.push((bases, after));
                    }
                }
                paths = longer;
            }
        }
        // The first edit that completes the stretch (substitutions first,
        // then deletions and insertions by length), else the one that walks
        // furthest.
        let best = tried.iter().position(|t| t.3).or_else(|| {
            let furthest = tried.iter().map(|t| t.2).max()?;
            tried.iter().position(|t| t.2 == furthest)
        })?;
        let (change, walked, reached, complete) = tried.swap_remove(best);
        found.push((pos, change));
        if complete {
            return Some((found, walked, reached))
        }
        pattern = walked;
        pos = reached;
    }
}

/// Calls the edits (SNPs and small indels) of `sequence` against the index
/// by walking its de Bruijn graph; sorted by position.
///
/// Where the contig's match in the index (its matching statistics), at least
/// as long as a random match is unlikely to be ([Index::random_match_threshold] for [MAX_ERROR_PROB]),
/// cannot be extended by the next base at p, every other base b is tried:
/// the match followed by b and then by the contig's bases, one step right in
/// the index at a time. If that walk takes the next _k_ − 1 bases, every
/// _k_-mer through p is in the index once p is b: a SNP. If no b gets that
/// far, the one that walks furthest is taken and the base where its walk
/// stops is substituted the same way, at most [MAX_EDITS] times. If
/// substitutions do not get the walk through, the same search is tried with
/// indels too: removing up to [MAX_DELETED] contig bases or adding up to
/// [MAX_INSERTED] index bases (every path of the index after the match);
/// a stretch whose walk does not get through has no edits. After an edit the
/// edited contig is walked on until its matching statistics are full _k_-mers
/// again, and a base not in the index there is repaired from that walk.
pub fn call_edits(
    sequence: &[u8],
    index: &Index,
) -> Vec<CalledEdit> {
    let k = index.k();
    let n = sequence.len();
    if n <= k {
        return Vec::new()
    }
    let min_match = index.random_match_threshold(MAX_ERROR_PROB);
    let ms = index.streaming().matching_statistics(sequence);

    let mut edits: Vec<Edit> = Vec::new();
    // The edited contig's match ending at i - 1, if a walk after an edit
    // stopped at i.
    let mut carried: Option<Pattern> = None;
    let mut i = 1;
    while i < n {
        let pattern = match carried.take() {
            Some(pattern) => pattern,
            None => {
                let (prev_len, len) = (ms[i - 1].0, ms[i].0);
                if prev_len < min_match || len >= (prev_len + 1).min(k) || base_code(sequence[i]).is_none() {
                    i += 1;
                    continue
                }
                (ms[i - 1].1.clone(), prev_len)
            },
        };

        let repaired = repair(index, sequence, pattern.clone(), i, false)
            .or_else(|| repair(index, sequence, pattern, i, true));
        match repaired {
            Some((found, mut pattern, mut pos)) => {
                for (at, change) in found {
                    let pos = at.try_into().expect("contig position does not fit in u32");
                    edits.push(match change {
                        TriedChange::Sub(b) => Edit { pos, change: Change::Sub { base: sequence[at], index_base: b } },
                        TriedChange::Del(d) => Edit { pos, change: Change::Del { bases: sequence[at..at + d].to_vec() } },
                        TriedChange::Ins(bases) => Edit { pos, change: Change::Ins { bases } },
                    });
                }
                // The original matching statistics still contain the edits here.
                while pos < n && ms[pos].0 < k {
                    match index.step_right(&pattern, sequence[pos]) {
                        Some(next) => {
                            pattern = next;
                            pos += 1;
                        },
                        None => break,
                    }
                }
                if pos < n && ms[pos].0 < k {
                    carried = Some(pattern);
                }
                i = pos;
            },
            None => i += 1,
        }
    }
    edits.into_iter().map(|edit| CalledEdit { key: key(sequence, &edit), edit }).collect()
}

/// `sequence` with `edits` (in its coordinates, increasing and not
/// overlapping) applied, and the edits at their positions in the result.
/// None if an edit does not match the sequence.
pub fn apply_edits(
    sequence: &[u8],
    edits: Vec<Edit>,
) -> Option<(Vec<u8>, Vec<Edit>)> {
    let mut out: Vec<u8> = Vec::with_capacity(sequence.len() + 16);
    let mut moved: Vec<Edit> = Vec::with_capacity(edits.len());
    let mut from = 0;
    for edit in edits {
        let pos = edit.pos as usize;
        if pos < from || pos > sequence.len() {
            return None
        }
        out.extend_from_slice(&sequence[from..pos]);
        let at = out.len() as u32;
        match &edit.change {
            Change::Sub { base, index_base } => {
                if sequence.get(pos) != Some(base) {
                    return None
                }
                out.push(*index_base);
                from = pos + 1;
            },
            Change::Del { bases } => {
                if sequence.get(pos..pos + bases.len()) != Some(bases.as_slice()) {
                    return None
                }
                from = pos + bases.len();
            },
            Change::Ins { bases } => {
                out.extend_from_slice(bases);
                from = pos;
            },
        }
        moved.push(Edit { pos: at, change: edit.change });
    }
    out.extend_from_slice(&sequence[from..]);
    Some((out, moved))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::BuildOpts;
    use crate::testing::random_sequence;

    #[test]
    fn edits_are_called_by_walking_the_index() {
        let reference = random_sequence(4000, 11);
        let k = 63;
        let index = Index::build(&[reference.clone()], &BuildOpts { k, ..Default::default() });

        // An isolated SNP, three within k of each other, one too close to the
        // end to be followed by k - 1 bases, an extra base in the contig (a
        // deletion), a missing one (an insertion), and a deleted base next to
        // a SNP (substitutions alone do not get through).
        let mut contig = reference.clone();
        let mut expected = Vec::new();
        for pos in [600, 1500, 1510, 1540, 3990] {
            let base = b"ACGT"[(base_code(reference[pos]).unwrap() as usize + 1) % 4];
            contig[pos] = base;
            if pos < 3990 {
                expected.push(Edit { pos: pos as u32, change: Change::Sub { base, index_base: reference[pos] } });
            }
        }
        let other = |pos: usize, by: usize| b"ACGT"[(base_code(reference[pos]).unwrap() as usize + by) % 4];
        let (extra, near) = (other(2499, 2), other(3300, 1));
        contig[3300] = near;
        contig.insert(3301, other(3301, 2));
        contig.remove(2900);
        contig.insert(2500, extra);
        let called: Vec<Edit> = call_edits(&contig, &index).into_iter().map(|called| called.edit).collect();
        let subs: Vec<Edit> = called.iter().filter(|e| matches!(e.change, Change::Sub { .. })).cloned().collect();
        assert_eq!(subs[..4], expected[..]);
        // The indels: the extra base and the missing one, and two edits for the
        // SNP next to an extra base (which two depends on the bases).
        assert!(called.iter().any(|e| matches!(e.change, Change::Del { .. })), "{:?}", called);
        assert!(called.iter().any(|e| matches!(e.change, Change::Ins { .. })), "{:?}", called);
        assert_eq!(called.len(), expected.len() + 4, "{:?}", called);
        // Editing gives the reference but for the SNP too close to the end.
        let (edited, _) = apply_edits(&contig, called).unwrap();
        let mut target = reference.clone();
        target[3990] = contig[contig.len() - 10];
        assert_eq!(edited, target);
        assert!(call_edits(&reference, &index).is_empty());
    }
}
