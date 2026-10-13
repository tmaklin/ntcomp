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

    /// Total number of recorded choices.
    pub fn n_choices(&self) -> usize {
        self.streams.iter().map(|x| x.n_choices as usize).sum()
    }

    /// Reader for the stream of `color`.
    pub fn reader(
        &self,
        color: u32,
    ) -> Option<ChoiceReader<'_>> {
        self.streams.get(color as usize).map(|stream| ChoiceReader { stream, bits: BitReader::new(&stream.bits), n_read: 0 })
    }
}

/// Choices of one colour.
#[derive(serde::Serialize, serde::Deserialize, Debug, Default, Clone, PartialEq, Eq)]
pub struct ChoiceStream {
    /// Number of choices in the stream.
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
    n_read: u32,
}

impl ChoiceReader<'_> {
    /// Reads the next choice among `n_options > 1` edges, None if the stream is exhausted.
    pub fn next(
        &mut self,
        n_options: u32,
    ) -> Option<u32> {
        if self.n_read == self.stream.n_choices {
            return None
        }
        let rank = self.bits.read_truncated_binary(n_options)?;
        self.n_read += 1;
        Some(rank)
    }

    /// True if every choice in the stream has been read.
    pub fn is_finished(&self) -> bool {
        self.n_read == self.stream.n_choices
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

