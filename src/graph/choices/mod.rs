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
//! Per-colour choice streams.
//!
//! Whenever a colour leaves a node with `n > 1` outgoing edges, the encoder
//! records the rank of the taken edge among all `n` edges with a truncated
//! binary code for `n` symbols, so a node with two outgoing edges costs one
//! bit. The decoder knows `n` at every step, so the codes need no separators.
//!
//! Each colour has its own stream so any contig can be decoded on its own.
//! A stream stores only its bits: the decoder stops at the sink, and then at
//! most the zero bits that pad the last byte are left.
//!
//! ## Compression
//!
//! The edges of every node are ordered by how many times the colours take
//! them (see [GraphBuilder](super::GraphBuilder)), so rank 0 is the most common
//! choice. Related colours also make the same choices through the same parts
//! of the graph. In the archive the choices of all colours are coded together
//! by a context-mixing model that predicts each colour's choices from those
//! of earlier colours ([model]); decoding needs the graph and
//! rebuilds these streams.
//!
use crate::codecs::bits::{BitReader, BitWriter};
pub mod model;

use super::ChainCsr;

/// Choice streams of every colour.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Choices {
    pub streams: Vec<ChoiceStream>,
}

impl Choices {
    /// Codes the streams of the colours walking `graph`, see [model].
    pub fn to_stored(
        &self,
        graph: &ChainCsr,
    ) -> Result<Vec<u8>, String> {
        model::encode(graph, self).map_err(|e| format!("choice streams: {}", e))
    }

    /// Decodes the streams of `n_colors` colours walking `graph`, see [Choices::to_stored].
    pub fn from_stored(
        stored: &[u8],
        graph: &ChainCsr,
        n_colors: usize,
    ) -> Result<Self, String> {
        model::decode(graph, n_colors, stored).map_err(|e| format!("choice streams: {}", e))
    }

    /// Total number of recorded choices.
    pub fn n_choices(&self) -> usize {
        self.streams.iter().map(|x| x.n_choices as usize).sum()
    }

    /// Reader for the stream of `color`.
    pub fn reader(
        &self,
        color: u32,
    ) -> Option<ChoiceReader<'_>> {
        self.streams.get(color as usize).map(|stream| ChoiceReader { stream, bits: BitReader::new(&stream.bits) })
    }
}

/// Choices of one colour.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct ChoiceStream {
    /// Number of choices in the stream, as written by [ChoiceWriter].
    pub n_choices: u32,

    /// Truncated binary codes, most significant bit first.
    pub bits: Vec<u8>,
}

#[derive(Default)]
pub struct ChoiceWriter {
    bits: BitWriter,
    n_choices: u32,
}

impl ChoiceWriter {
    pub fn new() -> Self {
        Self::default()
    }

    /// Records that edge `rank` of `n_options > 1` edges was taken.
    pub fn push(
        &mut self,
        rank: u32,
        n_options: u32,
    ) {
        assert!(n_options > 1 && rank < n_options, "invalid choice {} of {}", rank, n_options);
        self.bits.write_truncated_binary(rank, n_options);
        self.n_choices += 1;
    }

    pub fn finish(self) -> ChoiceStream {
        ChoiceStream { n_choices: self.n_choices, bits: self.bits.finish() }
    }
}

pub struct ChoiceReader<'a> {
    stream: &'a ChoiceStream,
    bits: BitReader<'a>,
}

impl ChoiceReader<'_> {
    /// Reads the next choice among `n_options > 1` edges, None if the stream is exhausted.
    pub fn next(
        &mut self,
        n_options: u32,
    ) -> Option<u32> {
        self.bits.read_truncated_binary(n_options)
    }

    /// True if what is left of the stream are fewer than 8 zero bits (the
    /// padding of the last byte). Meant for the check after the walk reached
    /// the sink: before that, choices whose codes are only zero bits and fit
    /// in the padding cannot be told from it.
    pub fn is_finished(&self) -> bool {
        let total = self.stream.bits.len() * 8;
        let read = self.bits.bits_read();
        if read > total || total - read >= 8 {
            return false
        }
        let left = total - read;
        left == 0 || self.stream.bits.last().is_some_and(|&b| b & ((1_u8 << left) - 1) == 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode(choices: &[Vec<(u32, u32)>]) -> Choices {
        let streams = choices.iter().map(|colour| {
            let mut writer = ChoiceWriter::new();
            for &(rank, n) in colour {
                writer.push(rank, n);
            }
            writer.finish()
        }).collect();
        Choices { streams }
    }

    #[test]
    fn choices_roundtrip() {
        // Every rank for 2..70 options, a large node and an empty colour.
        let all: Vec<(u32, u32)> = (2..70).flat_map(|n| (0..n).map(move |rank| (rank, n))).collect();
        let large = vec![(0, 100_000), (99_999, 100_000), (65_536, 70_000), (3, 4)];
        let choices = vec![all, vec![], large];
        let stored = encode(&choices);
        assert_eq!(stored.n_choices(), choices.iter().map(|c| c.len()).sum::<usize>());

        // Decode in reverse order to check that the colours are independent.
        for (color, expected) in choices.iter().enumerate().rev() {
            let mut reader = stored.reader(color as u32).unwrap();
            for (i, &(rank, n)) in expected.iter().enumerate() {
                // Not finished while choices are left.
                assert!(!reader.is_finished(), "colour {} finished after {} choices", color, i);
                assert_eq!(reader.next(n), Some(rank));
            }
            assert!(reader.is_finished());
        }
        assert!(stored.reader(3).is_none());
    }

    #[test]
    fn two_options_cost_one_bit() {
        let mut writer = ChoiceWriter::new();
        for i in 0..16 {
            writer.push(i % 2, 2);
        }
        assert_eq!(writer.finish().bits, vec![0b01010101, 0b01010101]);
    }
}
