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
//! Bit-level codes, written most significant bit first.
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

    /// Number of bits read so far.
    pub fn bits_read(&self) -> usize {
        self.bit
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_roundtrip() {
        let truncated: Vec<(u32, u32)> = (1..70).flat_map(|n| (0..n).map(move |v| (v, n))).collect();

        let mut writer = BitWriter::new();
        for &(v, n) in &truncated {
            writer.write_truncated_binary(v, n);
        }
        writer.write_bits(0b101, 3);
        let bytes = writer.finish();

        let mut reader = BitReader::new(&bytes);
        for &(v, n) in &truncated {
            assert_eq!(reader.read_truncated_binary(n), Some(v));
        }
        assert_eq!(reader.read_bits(3), Some(0b101));
        assert_eq!(reader.read_bits(8), None);
    }

    #[test]
    fn small_values_are_short() {
        let mut writer = BitWriter::new();
        writer.write_truncated_binary(1, 2); // 1
        writer.write_truncated_binary(2, 3); // 11
        writer.write_truncated_binary(0, 3); // 0
        assert_eq!(writer.finish(), vec![0b1110_0000]);
    }
}
