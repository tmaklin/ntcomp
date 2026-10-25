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
//! NLZ: a sequence as a dictionary of _k_-mers of an SBWT index.
//!
//! The parse ([parse]) covers a sequence, from its end to its start, with
//! the entries of its [Dictionary]: each a _k_-mer of the index (its colex
//! rank) and the number of bases it covers ([Entry]), the _k_-mer's last
//! bases or the _k_-mer left-extended through the de Bruijn graph of the
//! index for as long as the extension is unique. [Dictionary::decode]
//! recovers the sequence from the entries with the index. NLZ-OG ([crate::graph]) stores the dictionaries of many
//! sequences as coloured paths in a shared graph; [legacy] is the plain NLZ
//! file format, which stores them one after another.
//!
use std::ops::Range;

use crate::index::{Index, Match};
#[cfg(doc)]
use crate::walk::WalkIndex;
use crate::walk::IndexWalk;

pub mod legacy;

type E = Box<dyn std::error::Error>;

/// One entry of an NLZ dictionary: the _k_-mer at `colex_rank` of the index,
/// covering the `len` bases of the sequence before the next entry: the
/// _k_-mer's last `len` bases if `len` ≤ _k_, else the _k_-mer left-extended
/// to `len` bases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Entry {
    pub colex_rank: u32,
    pub len: u32,
    /// The _k_-mer is left-extended until the extension stops where
    /// [left_extend_kmer] stops (the index has no unique left extension), so
    /// decoding recovers `len` by extending the _k_-mer as far as possible;
    /// `len` may then be unknown (0).
    pub implicit: bool,
}

/// The NLZ dictionary of a sequence ([parse]): its entries from the end of
/// the sequence to its start, covering its `sequence_len` bases.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Dictionary {
    pub sequence_len: usize,
    pub entries: Vec<Entry>,
}

/// Left extends `kmer_start` while the extension is unique in the SBWT, by at
/// most `max_extension_len` characters.
///
/// If `ref_nts` is given, also checks that the extension matches `ref_nts`.
///
/// Returns the extended k-mer and whether the extension stopped because the
/// SBWT has no unique left extension, i.e. where [left_extend_kmer] stops as
/// well. False if it stopped because of a mismatch with `ref_nts` or because
/// `max_extension_len` was reached.
///
/// Four searches per character: the reference the walks of
/// [IndexWalk::left_extend] (one step per
/// character, from a known row) are checked against, and the left extension
/// of the plain NLZ format ([legacy]).
pub fn left_extend_kmer(
    kmer_start: &[u8],
    ref_nts: Option<&[u8]>,
    sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    max_extension_len: usize,
) -> (Vec<u8>, bool) {
    assert!(!kmer_start.is_empty());
    if let Some(ref_nts) = ref_nts {
        assert!(!ref_nts.is_empty());
    }

    let mut left_extension_len = 0;
    let mut kmer = kmer_start.to_vec();
    let mut stopped_in_sbwt = false;
    while left_extension_len < max_extension_len {
        let new_kmers: Vec<(Vec<u8>, Range<usize>)> = sbwt.alphabet().iter().filter_map(|c| {
            let new_kmer: Vec<u8> = [&[*c], &kmer[0..(kmer.len() - (left_extension_len + 1))]].concat();
            let res = sbwt.search(&new_kmer);
            if res.as_ref().is_some() {
                Some((new_kmer, res.unwrap()))
            } else {
                None
            }
        }).collect();
        if !new_kmers.is_empty() {
            let unique = new_kmers.len() == 1 && new_kmers[0].1.end - new_kmers[0].1.start == 1;

            let seq_matches = if let Some(ref_nts) = ref_nts {
                new_kmers[0].0[0] == ref_nts[ref_nts.len() - kmer.len() - 1]
            } else {
                true
            };

            if seq_matches && unique {
                kmer = [&[new_kmers[0].0[0]], kmer.as_slice()].concat();
            } else {
                stopped_in_sbwt = !unique;
                break;
            }
        } else {
            stopped_in_sbwt = true;
            break;
        }
        left_extension_len += 1;
    }
    (kmer, stopped_in_sbwt)
}

/// The NLZ parse of a nucleotide sequence against an SBWT index: its
/// [Dictionary], from the end of the sequence to the start.
pub fn parse(
    nucleotides: &[u8],
    index: &Index,
) -> Result<Dictionary, E> {
    parse_with(nucleotides, index.matching_statistics(nucleotides), index)
}

/// Same as [parse] with the matching statistics `ms` of `nucleotides`
/// ([Index::matching_statistics]).
pub fn parse_with(
    nucleotides: &[u8],
    ms: Vec<Match>,
    index: &Index,
) -> Result<Dictionary, E> {
    let n = nucleotides.len();
    let mut implicit: Vec<bool> = Vec::new();
    let matches: Vec<(usize, usize)> = {
        let k = index.k();
        let mut res = ms;

        let mut i = n;
        let mut kept: Vec<usize> = Vec::new();
        while i > 0 {
            if res[i - 1].len as usize == k && i > k + 1 {
                // Left extend as [left_extend_kmer] with the sequence as the
                // reference, through the rows of the matching statistics: the
                // k-mer ending one base to the left is the predecessor the
                // sequence takes, and its group holds all the predecessors.
                let max_extension_len = i - k - 1;
                let mut row = res[i - 1].start as usize;
                let mut extension_len = 0;
                let mut stopped_in_sbwt = false;
                while extension_len < max_extension_len {
                    let left = i - 2 - extension_len;
                    if res[left].len as usize == k {
                        let predecessor = res[left].start as usize;
                        if !index.only_in_group(predecessor) {
                            stopped_in_sbwt = true;
                            break
                        }
                        row = predecessor;
                        extension_len += 1;
                    } else {
                        // The sequence's predecessor is not in the index: the
                        // extension stops in the SBWT unless the k-mer has exactly
                        // one (real) predecessor, a different one.
                        stopped_in_sbwt = !index.sbwt().inverse_lf_step(row).is_some_and(|predecessor| {
                            index.only_in_group(predecessor) && index.is_real_kmer(predecessor)
                        });
                        break
                    }
                }
                let phrase_len = k + extension_len;

                // Store the whole left extension; the next entry ends where it starts.
                res[i - 1].len = phrase_len as u32;
                kept.push(i - 1);
                implicit.push(stopped_in_sbwt);
                i -= phrase_len - 1;
            } else {
                kept.push(i - 1);
                implicit.push(false);
                if i > res[i - 1].len as usize {
                    i -= res[i - 1].len as usize - 1;
                } else {
                    break;
                }
            }

            if i > 0 {
                i -= 1;
            } else {
                break;
            }
        }
        kept.iter().map(|x| (res[*x].len as usize, res[*x].start as usize)).collect()
    };

    let dictionary_bases = matches.iter().map(|x| x.0).sum::<usize>();

    assert_eq!(dictionary_bases, nucleotides.len());
    assert_eq!(matches.len(), implicit.len());

    let entries = matches.into_iter().zip(implicit).map(|((len, colex_rank), implicit)| {
        let colex_rank = u32::try_from(colex_rank).map_err(|_| "colex rank does not fit in u32")?;
        let len = u32::try_from(len).map_err(|_| "dictionary entry longer than 2^32 bases")?;
        Ok(Entry { colex_rank, len, implicit })
    }).collect::<Result<Vec<Entry>, &str>>()?;
    Ok(Dictionary { sequence_len: n, entries })
}

impl Dictionary {
    /// Decodes the sequence.
    ///
    /// The entries are processed in their own order, from the end of the
    /// sequence to the start, so that the end position of every entry is
    /// known. An [implicit](Entry::implicit) entry's _k_-mer is extended as
    /// far as the encoder could have extended it at that position; any other
    /// is extended to its `len` if `len` > _k_, and is the _k_-mer's last
    /// `len` bases otherwise.
    ///
    /// The _k_-mers are extended by walking the de Bruijn graph of the
    /// `index` ([IndexWalk::left_extend]): an [Index], or a [WalkIndex]
    /// without the LCS array (the LCS values are then computed from the
    /// SBWT).
    ///
    /// Fails if the entries do not cover exactly `sequence_len` bases.
    pub fn decode<I: IndexWalk>(
        &self,
        index: &I,
    ) -> Result<Vec<u8>, DecodeError> {
        let sequence_length = self.sequence_len;
        let mut sequence: Vec<u8> = vec![0; sequence_length];
        let mut end = sequence_length;
        let (sbwt, k) = (index.sbwt(), index.k());
        // Left extension of the k-mer at `colex_rank` by walking the SBWT.
        let extend = |colex_rank: u32, max_extension_len: usize, implicit: bool| index.left_extend(colex_rank as usize, max_extension_len, implicit).0;
        for entry in &self.entries {
            let (kmer, len) = if entry.implicit {
                // The encoder extends at most end - k - 1 characters.
                let new_kmer = extend(entry.colex_rank, end.saturating_sub(k + 1), true);
                let len = new_kmer.len();
                (new_kmer, len)
            } else if entry.len as usize > k {
                let new_kmer = extend(entry.colex_rank, entry.len as usize - k, false);
                assert_eq!(new_kmer.len(), entry.len as usize);
                (new_kmer, entry.len as usize)
            } else {
                (sbwt.access_kmer(entry.colex_rank as usize), entry.len as usize)
            };
            if len > end {
                return Err(DecodeError(format!("dictionary covers more than {} bases", sequence_length)))
            }
            sequence[(end - len)..end].copy_from_slice(&kmer[(kmer.len() - len)..]);
            end -= len;
        }
        if end != 0 {
            return Err(DecodeError(format!("dictionary covers {} of {} bases", sequence_length - end, sequence_length)))
        }

        Ok(sequence)
    }
}

/// Error from [Dictionary::decode]: the dictionary does not cover the
/// sequence.
#[derive(Debug)]
pub struct DecodeError(pub String);

impl std::fmt::Display for DecodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "NLZ dictionary: {}", self.0)
    }
}

impl std::error::Error for DecodeError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::related_index;

    #[test]
    fn left_extensions_with_and_without_lcs_agree() {
        // Related sequences with substitutions, so the de Bruijn graph branches.
        let (seqs, index) = related_index(15);
        let (sbwt, k) = (index.sbwt(), index.k());

        let mut n_extended = 0;
        for seq in &seqs {
            for end in k..=seq.len() {
                let kmer = &seq[end - k..end];
                let max = end - k;
                let searched = left_extend_kmer(kmer, None, sbwt, max);
                let plain = searched.0.clone();
                // Walking the SBWT from the k-mer's row, with the LCS array or the
                // LCS values from the SBWT: checked, and by a known length.
                let row = sbwt.search(kmer).unwrap().start;
                let (with_array, without) = (&index, crate::walk::WalkIndex { sbwt, lcs: None });
                assert_eq!(with_array.left_extend(row, max, true), searched);
                assert_eq!(without.left_extend(row, max, true), searched);
                assert_eq!(with_array.left_extend(row, plain.len() - k, false).0, plain);
                n_extended += (plain.len() > k) as usize;
            }

            // The encoder's extensions from the matching statistics are those of
            // left_extend_kmer with the sequence as the reference.
            let dictionary = parse(seq, &index).unwrap();
            assert_eq!(dictionary.sequence_len, seq.len());
            let mut end = seq.len();
            for entry in &dictionary.entries {
                let len = entry.len as usize;
                if len >= k && end > k + 1 {
                    let extended = left_extend_kmer(&seq[end - k..end], Some(&seq[..end]), sbwt, end - k - 1);
                    assert_eq!(extended, (seq[end - len..end].to_vec(), entry.implicit));
                }
                end -= len;
            }

            // Both decoders recover the sequence from its dictionary, also with
            // the implicit lengths unknown (as the graph stores them).
            let stored = Dictionary {
                entries: dictionary.entries.iter().map(|&e| Entry { len: if e.implicit { 0 } else { e.len }, ..e }).collect(),
                ..dictionary.clone()
            };
            for dictionary in [&dictionary, &stored] {
                assert_eq!(&dictionary.decode(&crate::walk::WalkIndex { sbwt, lcs: None }).unwrap(), seq);
                assert_eq!(&dictionary.decode(&index).unwrap(), seq);
            }
        }
        assert!(n_extended > 0, "no k-mer was left extended");
    }
}
