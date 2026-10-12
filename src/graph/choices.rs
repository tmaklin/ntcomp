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

#[derive(serde::Serialize, serde::Deserialize, Debug, Default)]
pub struct Choices {
    pub streams: Vec<ChoiceStream>,
}

impl Choices {
    /// Total number of recorded choices.
    pub fn n_choices(&self) -> usize {
        self.streams.iter().map(|x| x.n_choices as usize).sum()
    }

    /// Reader for the stream of `color`.
    pub fn reader(
        &self,
        color: u32,
    ) -> Option<ChoiceReader<'_>> {
        self.streams.get(color as usize).map(|stream| ChoiceReader { stream, bit: 0, n_read: 0 })
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

/// Number of bits `k` and threshold `u` of the truncated binary code for `n` symbols.
///
/// Ranks below `u` are written with `k` bits, the rest as `rank + u` with `k + 1` bits.
#[inline]
fn truncated_binary(
    n_options: u32,
) -> (u32, u32) {
    let k = 31 - n_options.leading_zeros();
    let u = ((2_u64 << k) - n_options as u64) as u32;
    (k, u)
}

#[derive(Default)]
pub struct ChoiceWriter {
    stream: ChoiceStream,
    n_bits: usize,
}

impl ChoiceWriter {
    pub fn new() -> Self {
        Self::default()
    }

    fn write_bits(
        &mut self,
        value: u32,
        n_bits: u32,
    ) {
        for i in (0..n_bits).rev() {
            if self.n_bits % 8 == 0 {
                self.stream.bits.push(0);
            }
            if (value >> i) & 1 == 1 {
                *self.stream.bits.last_mut().unwrap() |= 0x80 >> (self.n_bits % 8);
            }
            self.n_bits += 1;
        }
    }

    /// Records that edge `rank` of `n_options > 1` edges was taken.
    pub fn push(
        &mut self,
        rank: u32,
        n_options: u32,
    ) {
        assert!(n_options > 1 && rank < n_options, "invalid choice {} of {}", rank, n_options);
        let (k, u) = truncated_binary(n_options);
        if rank < u {
            self.write_bits(rank, k);
        } else {
            self.write_bits(rank + u, k + 1);
        }
        self.stream.n_choices += 1;
    }

    pub fn finish(self) -> ChoiceStream {
        self.stream
    }
}

pub struct ChoiceReader<'a> {
    stream: &'a ChoiceStream,
    bit: usize,
    n_read: u32,
}

impl ChoiceReader<'_> {
    fn read_bits(
        &mut self,
        n_bits: u32,
    ) -> Option<u32> {
        if self.bit + n_bits as usize > self.stream.bits.len() * 8 {
            return None
        }
        let mut value = 0;
        for _ in 0..n_bits {
            let bit = (self.stream.bits[self.bit / 8] >> (7 - self.bit % 8)) & 1;
            value = (value << 1) | bit as u32;
            self.bit += 1;
        }
        Some(value)
    }

    /// Reads the next choice among `n_options > 1` edges, None if the stream is exhausted.
    pub fn next(
        &mut self,
        n_options: u32,
    ) -> Option<u32> {
        if self.n_read == self.stream.n_choices {
            return None
        }
        let (k, u) = truncated_binary(n_options);
        let value = self.read_bits(k)?;
        let rank = if value < u {
            value
        } else {
            ((value << 1) | self.read_bits(1)?) - u
        };
        self.n_read += 1;
        Some(rank)
    }

    /// True if every choice in the stream has been read.
    pub fn is_finished(&self) -> bool {
        self.n_read == self.stream.n_choices
    }
}
