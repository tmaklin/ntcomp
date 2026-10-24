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
//! entries that are each a _k_-mer of the index (its colex rank) and a suffix
//! length ([SuffixLen]): the entry's last bases, or the _k_-mer left-extended
//! through the de Bruijn graph of the index for as long as the extension is
//! unique. [decode_dictionary] recovers the sequence from the entries with
//! the index. NLZ-OG ([crate::graph]) stores the dictionaries of many
//! sequences as coloured paths in a shared graph; [legacy] is the plain NLZ
//! file format, which stores them one after another.
//!
use std::ops::Range;

use sbwt::LcsArray;
use sbwt::StreamingIndex;
use sbwt::sbwt_index_variant::SbwtIndexVariant;

use crate::walk::{IndexWalk, LcsIndex, WalkIndex};

pub mod legacy;

type E = Box<dyn std::error::Error>;

/// [Implicit](SuffixLen::Implicit) suffix lengths are stored as a 0; a real
/// suffix length can never be 0.
pub const IMPLICIT_SUFFIX_LEN: u32 = 0;

/// Suffix length and type for an edge
///
/// Classifies the suffix length for a dictionary entry by specifying whether
/// its implicit (for a left extended _k_-mer but can be inferred), explicit
/// (left extended and needs to be declared), or plain (subsequence of a plan
/// _k_-mer). Stores the value for explicit and plain types.
///
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SuffixLen {
    /// The dictionary entry is a left extended _k_-mer whose extension stopped
    /// where [left_extend_kmer] stops, so [decode_dictionary] recovers the
    /// length by extending the _k_-mer as far as possible.
    Implicit,

    /// Entry is a left extended _k_-mer, whose extension stops after an
    /// explicit number of extensions, so [decode_dictionary] needs to know the
    /// length.
    Explicit { suffix_len: u32 },

    /// Entry is a plain _k_-mer whose suffix length is less or equal to the
    /// value of k.
    Plain { suffix_len: u32 },
}

impl SuffixLen {
    /// The suffix length `len` as stored on the graph's edges
    /// ([IMPLICIT_SUFFIX_LEN] for an implicit one), for _k_-mer size `k`.
    pub fn from_stored(
        len: u32,
        k: u32,
    ) -> SuffixLen {
        if len == IMPLICIT_SUFFIX_LEN {
            SuffixLen::Implicit
        } else if len > k {
            SuffixLen::Explicit { suffix_len: len }
        } else {
            SuffixLen::Plain { suffix_len: len }
        }
    }
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
/// dictionary, from the end of the sequence to the start.
///
/// Returns the dictionary and, for every dictionary entry, a boolean denoting
/// whether its suffix length is *implicit*: the entry is a left extended k-mer
/// whose extension stopped where [left_extend_kmer] stops as well, so the
/// decoder can recover the length by extending the k-mer as far as possible.
pub fn parse(
    nucleotides: &[u8],
    sbwt: &SbwtIndexVariant,
    lcs: &LcsArray,
) -> Result<(Vec<(usize, Range<usize>)>, Vec<bool>), E> {
    let n = nucleotides.len();
    let mut implicit: Vec<bool> = Vec::new();
    let dictionary: Vec<(usize, Range<usize>)> = match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            let index = StreamingIndex::new(sbwt, lcs);
            let mut res = index.matching_statistics(nucleotides);
            let walk = LcsIndex { sbwt, lcs };

            let mut i = n;
            let mut kept: Vec<usize> = Vec::new();
            while i > 0 {
                if res[i - 1].0 == k && i > k + 1 {
                    // Left extend as [left_extend_kmer] with the sequence as the
                    // reference, through the rows of the matching statistics: the
                    // k-mer ending one base to the left is the predecessor the
                    // sequence takes, and its group holds all the predecessors.
                    let max_extension_len = i - k - 1;
                    let mut row = res[i - 1].1.start;
                    let mut extension_len = 0;
                    let mut stopped_in_sbwt = false;
                    while extension_len < max_extension_len {
                        let left = i - 2 - extension_len;
                        if res[left].0 == k {
                            let predecessor = res[left].1.start;
                            if !walk.only_in_group(predecessor) {
                                stopped_in_sbwt = true;
                                break
                            }
                            row = predecessor;
                            extension_len += 1;
                        } else {
                            // The sequence's predecessor is not in the index: the
                            // extension stops in the SBWT unless the k-mer has exactly
                            // one (real) predecessor, a different one.
                            stopped_in_sbwt = !sbwt.inverse_lf_step(row).is_some_and(|predecessor| {
                                walk.only_in_group(predecessor) && walk.is_real_kmer(predecessor)
                            });
                            break
                        }
                    }
                    let phrase_len = k + extension_len;

                    // Store the whole left extension; the next entry ends where it starts.
                    res[i - 1].0 = phrase_len;
                    kept.push(i - 1);
                    implicit.push(stopped_in_sbwt);
                    i -= phrase_len - 1;
                } else {
                    kept.push(i - 1);
                    implicit.push(false);
                    if i > res[i - 1].0 {
                        i -= res[i - 1].0 - 1;
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
            kept.iter().map(|x| res[*x].clone()).collect()
        },
    };

    let dictionary_bases = dictionary.iter().map(|x| x.0).sum::<usize>();

    assert_eq!(dictionary_bases, nucleotides.len());
    assert_eq!(dictionary.len(), implicit.len());

    Ok((dictionary, implicit))
}

/// Decodes the nucleotide sequence of length `sequence_length` from its NLZ
/// dictionary (colex rank, suffix length): a path of the graph
/// ([GraphWalk::walk_path](crate::graph::GraphWalk::walk_path), with
/// [SuffixLen::from_stored]).
///
/// The dictionary is processed in its own order, from the end of the sequence
/// to the start, so that the end position of every entry is known.
/// [Implicit](SuffixLen::Implicit) suffix lengths are recovered by extending
/// the k-mer as far as the encoder could have extended it at that position and
/// explicit and plain by taking the suffix length stored inside [SuffixLen].
///
/// The k-mers are extended by walking the de Bruijn graph of the SBWT
/// ([IndexWalk::left_extend](crate::walk::IndexWalk::left_extend)), with
/// `lcs` if given and with LCS values computed from the SBWT otherwise.
///
/// Fails if the dictionary does not cover exactly `sequence_length` bases.
pub fn decode_dictionary(
    dictionary: Vec<(u32, SuffixLen)>,
    sequence_length: usize,
    sbwt: &SbwtIndexVariant,
    lcs: Option<&LcsArray>,
) -> Result<Vec<u8>, DecodeError> {
    let mut sequence: Vec<u8> = vec![0; sequence_length];
    let mut end = sequence_length;
    match sbwt {
        SbwtIndexVariant::SubsetMatrix(sbwt) => {
            let k = sbwt.k();
            // Left extension of the k-mer at `colex_rank` by walking the SBWT.
            let walk = WalkIndex { sbwt, lcs };
            let extend = |colex_rank: u32, max_extension_len: usize, implicit: bool| walk.left_extend(colex_rank as usize, max_extension_len, implicit).0;
            for (colex_rank, suffix_len) in dictionary {
                let (kmer, len) = match suffix_len {
                    SuffixLen::Implicit => {
                        // The encoder extends at most end - k - 1 characters.
                        let new_kmer = extend(colex_rank, end.saturating_sub(k + 1), true);
                        let len = new_kmer.len();
                        (new_kmer, len)
                    },
                    SuffixLen::Explicit { suffix_len: len } => {
                        let new_kmer = extend(colex_rank, (len - k as u32) as usize, false);
                        assert_eq!(new_kmer.len(), len as usize);
                        (new_kmer, len as usize)
                    },
                    SuffixLen::Plain { suffix_len: len } => {
                        (sbwt.access_kmer(colex_rank as usize), len as usize)
                    },
                };
                if len > end {
                    return Err(DecodeError(format!("dictionary covers more than {} bases", sequence_length)))
                }
                sequence[(end - len)..end].copy_from_slice(&kmer[(kmer.len() - len)..]);
                end -= len;
            }
        },
    };
    if end != 0 {
        return Err(DecodeError(format!("dictionary covers {} of {} bases", sequence_length - end, sequence_length)))
    }

    Ok(sequence)
}

/// Error from [decode_dictionary]: the dictionary does not cover the sequence.
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
        let (seqs, variant, lcs) = related_index(15);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = &variant;
        let k = sbwt.k();

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
                let (with_array, without) = (LcsIndex { sbwt, lcs: &lcs }, WalkIndex { sbwt, lcs: None });
                assert_eq!(with_array.left_extend(row, max, true), searched);
                assert_eq!(without.left_extend(row, max, true), searched);
                assert_eq!(with_array.left_extend(row, plain.len() - k, false).0, plain);
                n_extended += (plain.len() > k) as usize;
            }

            // The encoder's extensions from the matching statistics are those of
            // left_extend_kmer with the sequence as the reference.
            let (dictionary, implicit) = parse(seq, &variant, &lcs).unwrap();
            let mut end = seq.len();
            for ((len, _), &imp) in dictionary.iter().zip(&implicit) {
                if *len >= k && end > k + 1 {
                    let extended = left_extend_kmer(&seq[end - k..end], Some(&seq[..end]), sbwt, end - k - 1);
                    assert_eq!(extended, (seq[end - len..end].to_vec(), imp));
                }
                end -= len;
            }

            // Both decoders recover the sequence from its encoding.
            let path = |_| dictionary.iter().zip(&implicit).map(|((len, colex), &imp)| {
                let len = *len as u32;
                let suffix_len = if imp {
                    SuffixLen::Implicit
                } else if len > k as u32 {
                    SuffixLen::Explicit { suffix_len: len }
                } else {
                    SuffixLen::Plain { suffix_len: len }
                };
                (colex.start as u32, suffix_len)
            }).collect::<Vec<_>>();
            assert_eq!(&decode_dictionary(path(()), seq.len(), &variant, None).unwrap(), seq);
            assert_eq!(&decode_dictionary(path(()), seq.len(), &variant, Some(&lcs)).unwrap(), seq);
        }
        assert!(n_extended > 0, "no k-mer was left extended");
    }
}
