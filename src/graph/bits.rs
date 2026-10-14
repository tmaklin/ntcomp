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
//! Bit-level codes used in the stored graph, written most significant bit first.
//!

/// Number of bits `k` and threshold `u` of the truncated binary code for `n` symbols.
///
/// Values below `u` are written with `k` bits, the rest as `value + u` with `k + 1` bits.
fn truncated_binary(
    n_symbols: u32,
) -> (u32, u32) {
    let k = 31 - n_symbols.leading_zeros();
    let u = ((2_u64 << k) - n_symbols as u64) as u32;
    (k, u)
}

#[derive(Default)]
pub struct BitWriter {
    bytes: Vec<u8>,
    n_bits: usize,
}

impl BitWriter {
    pub fn new() -> Self {
        Self::default()
    }

    /// Writes the lowest `n_bits` bits of `value`.
    pub fn write_bits(
        &mut self,
        value: u64,
        n_bits: u32,
    ) {
        for i in (0..n_bits).rev() {
            if self.n_bits.is_multiple_of(8) {
                self.bytes.push(0);
            }
            if (value >> i) & 1 == 1 {
                *self.bytes.last_mut().unwrap() |= 0x80 >> (self.n_bits % 8);
            }
            self.n_bits += 1;
        }
    }

    /// Writes `value > 0` with the Elias gamma code.
    pub fn write_gamma(
        &mut self,
        value: u64,
    ) {
        assert!(value > 0, "Elias gamma code needs a positive value");
        let n = 63 - value.leading_zeros();
        self.write_bits(0, n);
        self.write_bits(value, n + 1);
    }

    /// Writes `value < n_symbols` with the truncated binary code for `n_symbols`.
    pub fn write_truncated_binary(
        &mut self,
        value: u32,
        n_symbols: u32,
    ) {
        assert!(value < n_symbols, "{} does not fit in {} symbols", value, n_symbols);
        let (k, u) = truncated_binary(n_symbols);
        if value < u {
            self.write_bits(value as u64, k);
        } else {
            self.write_bits((value + u) as u64, k + 1);
        }
    }

    pub fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

pub struct BitReader<'a> {
    bytes: &'a [u8],
    bit: usize,
}

impl<'a> BitReader<'a> {
    pub fn new(
        bytes: &'a [u8],
    ) -> Self {
        BitReader { bytes, bit: 0 }
    }

    /// Reads `n_bits <= 64` bits, None past the end of the input.
    pub fn read_bits(
        &mut self,
        n_bits: u32,
    ) -> Option<u64> {
        if n_bits > 64 || self.bit + n_bits as usize > self.bytes.len() * 8 {
            return None
        }
        let mut value = 0;
        for _ in 0..n_bits {
            let bit = (self.bytes[self.bit / 8] >> (7 - self.bit % 8)) & 1;
            value = (value << 1) | bit as u64;
            self.bit += 1;
        }
        Some(value)
    }

    pub fn read_gamma(
        &mut self,
    ) -> Option<u64> {
        let mut n = 0;
        while self.read_bits(1)? == 0 {
            n += 1;
            if n > 63 {
                return None
            }
        }
        Some((1 << n) | self.read_bits(n)?)
    }

    pub fn read_truncated_binary(
        &mut self,
        n_symbols: u32,
    ) -> Option<u32> {
        let (k, u) = truncated_binary(n_symbols);
        let value = self.read_bits(k)? as u32;
        if value < u {
            Some(value)
        } else {
            Some((((value as u64) << 1 | self.read_bits(1)?) - u as u64) as u32)
        }
    }
}

