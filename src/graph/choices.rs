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
//! of the graph. The streams of all colours are therefore stored together as
//! one LZMA stream ([LzmaData]), which is decompressed when the choices
//! are loaded. Decoding a colour only reads its own stream.
//!
use super::bits::{BitReader, BitWriter};
use super::lzma_wrapper::{self, LzmaData};

/// Choice streams of every colour.
///
/// Serialized as [LzmaData] holding the postcard serialized streams.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Choices {
    pub streams: Vec<ChoiceStream>,
}

impl Choices {
    /// Compresses the streams.
    pub fn to_stored(
        &self,
    ) -> Result<LzmaData, String> {
        let streams = postcard::to_allocvec(&self.streams).map_err(|e| format!("choice streams: {}", e))?;
        lzma_wrapper::compress(&streams).map_err(|e| format!("choice streams: {}", e))
    }

    /// Decompresses the streams.
    pub fn from_stored(
        stored: LzmaData,
    ) -> Result<Self, String> {
        let streams = lzma_wrapper::decompress(&stored).map_err(|e| format!("choice streams: {}", e))?;
        let streams = postcard::from_bytes(&streams).map_err(|e| format!("choice streams: {}", e))?;
        Ok(Choices { streams })
    }

    /// Total number of recorded choices (0 for decoded choices, see
    /// [ChoiceStream::n_choices]).
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
#[derive(serde::Serialize, serde::Deserialize, Debug, Default, Clone, PartialEq, Eq)]
pub struct ChoiceStream {
    /// Number of choices in the stream, as written by [ChoiceWriter]; not
    /// stored (0 in decoded streams).
    #[serde(skip)]
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

impl serde::Serialize for Choices {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.to_stored().map_err(serde::ser::Error::custom)?.serialize(serializer)
    }
}

impl<'de> serde::Deserialize<'de> for Choices {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let stored = LzmaData::deserialize(deserializer)?;
        Choices::from_stored(stored).map_err(serde::de::Error::custom)
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

        let bytes = postcard::to_allocvec(&stored).unwrap();
        let decoded: Choices = postcard::from_bytes(&bytes).unwrap();
        // The bits are stored, the counts are not.
        assert_eq!(decoded.streams.iter().map(|s| &s.bits).collect::<Vec<_>>(), stored.streams.iter().map(|s| &s.bits).collect::<Vec<_>>());
        assert!(decoded.streams.iter().all(|s| s.n_choices == 0));

        // Decode in reverse order to check that the colours are independent.
        for (color, expected) in choices.iter().enumerate().rev() {
            let mut reader = decoded.reader(color as u32).unwrap();
            for (i, &(rank, n)) in expected.iter().enumerate() {
                // Not finished while choices are left.
                assert!(!reader.is_finished(), "colour {} finished after {} choices", color, i);
                assert_eq!(reader.next(n), Some(rank));
            }
            assert!(reader.is_finished());
        }
        assert!(decoded.reader(3).is_none());
    }

    #[test]
    fn two_options_cost_one_bit() {
        let mut writer = ChoiceWriter::new();
        for i in 0..16 {
            writer.push(i % 2, 2);
        }
        assert_eq!(writer.finish().bits, vec![0b01010101, 0b01010101]);
    }

    #[test]
    fn repeated_choices_compress() {
        // Many colours making the same, mostly rank 0, choices.
        let colour: Vec<(u32, u32)> = (0..2000).map(|i| ((i % 13 == 0) as u32, 2 + (i % 3))).collect();
        let choices = encode(&vec![colour; 50]);
        let uncompressed = postcard::to_allocvec(&choices.streams).unwrap().len();
        let stored = postcard::to_allocvec(&choices).unwrap().len();
        assert!(stored * 20 < uncompressed, "{} vs {}", stored, uncompressed);
    }

    #[test]
    fn empty_and_corrupt_choices() {
        let empty = Choices::default();
        let decoded: Choices = postcard::from_bytes(&postcard::to_allocvec(&empty).unwrap()).unwrap();
        assert_eq!(decoded, empty);

        let colour: Vec<(u32, u32)> = (0..2000).map(|i| ((i * 7 % 5) as u32, 5)).collect();
        let mut stored = encode(&[colour]).to_stored().unwrap();
        stored.bytes.truncate(stored.bytes.len() / 2);
        assert!(Choices::from_stored(stored).is_err());
    }
}
