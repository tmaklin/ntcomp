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
//! Binary arithmetic coding with adaptive bit models, shared by the index
//! tig coder ([super::tigs]) and the archive format ([crate::graph]).
//!
//! [Encoder] and [Decoder] are a carryless 32-bit binary arithmetic coder (as
//! in lpaq) driven by 16-bit probabilities of a 1. Both implement [Coder], so
//! one function can both encode and decode a value: it passes the true bits
//! (ignored when decoding) and gets back the coded ones.
//!
//! An adaptive bit model is a [u32] counter: the probability of a 1 in the
//! high 22 bits and a count in the low 10. [update] moves the probability
//! towards each coded bit by 1 / (count + 1.5), so a model adapts fast at
//! first and settles as the count grows (up to [LIMIT]). Everything is integer
//! arithmetic, so encoder and decoder agree on every platform.

/// A fresh adaptive bit model: probability ½, count 0.
pub const C0: u32 = (1 << 21) << 10;

/// Largest count of an adaptive bit model (its slowest adaptation rate).
pub const LIMIT: u32 = 255;

/// 65536 / (n + 1.5), rounded down, for counts n = 0..1024.
const DIV: [u32; 1024] = {
    let mut t = [0u32; 1024];
    let mut n = 0;
    while n < 1024 { t[n] = 131_072 / (2 * n as u32 + 3); n += 1; }
    t
};

/// Updates an adaptive bit model with a coded `bit`.
#[inline]
pub fn update(
    counter: &mut u32,
    bit: u32,
) {
    let n = *counter & 1023;
    let p = (*counter >> 10) as i64;
    let target = if bit == 1 { (1i64 << 22) - 1 } else { 0 };
    let p = p + (((target - p) * DIV[n as usize] as i64) >> 16);
    *counter = ((p as u32) << 10) | (n + 1).min(LIMIT);
}

/// The 16-bit probability of a 1 of an adaptive bit model, within 1..=65535.
#[inline]
pub fn p1(
    counter: u32,
) -> u32 {
    (counter >> 16).clamp(1, 65535)
}

/// Arithmetic encoder.
pub struct Encoder { x1: u32, x2: u32, out: Vec<u8> }

impl Default for Encoder {
    fn default() -> Self { Self::new() }
}

impl Encoder {
    pub fn new() -> Encoder { Encoder { x1: 0, x2: u32::MAX, out: Vec::new() } }

    /// Codes `bit` with probability `p1` (16-bit, 1..=65535) of a 1.
    pub fn bit(&mut self, bit: u32, p1: u32) {
        let xmid = self.x1 + (((self.x2 - self.x1) as u64 * p1 as u64) >> 16) as u32;
        if bit == 1 { self.x2 = xmid } else { self.x1 = xmid + 1 }
        while (self.x1 ^ self.x2) & 0xFF00_0000 == 0 { self.out.push((self.x2 >> 24) as u8); self.x1 <<= 8; self.x2 = (self.x2 << 8) | 255; }
    }

    /// Bytes written so far.
    pub fn len(&self) -> usize { self.out.len() }

    /// Whether nothing has been written yet.
    pub fn is_empty(&self) -> bool { self.out.is_empty() }

    /// Flushes the coder and returns the coded bytes.
    pub fn finish(mut self) -> Vec<u8> { for _ in 0..4 { self.out.push((self.x1 >> 24) as u8); self.x1 <<= 8; } self.out }
}

/// Arithmetic decoder of bytes written by [Encoder]. Reading past the end
/// yields zero bytes (corrupt input decodes to garbage, not a panic) and is
/// reported by [Decoder::overrun].
pub struct Decoder<'a> { x1: u32, x2: u32, x: u32, inp: &'a [u8], pos: usize }

impl<'a> Decoder<'a> {
    pub fn new(inp: &'a [u8]) -> Decoder<'a> {
        let mut d = Decoder { x1: 0, x2: u32::MAX, x: 0, inp, pos: 0 };
        for _ in 0..4 { d.x = (d.x << 8) | d.byte(); }
        d
    }

    fn byte(&mut self) -> u32 { let b = self.inp.get(self.pos).copied().unwrap_or(0); self.pos += 1; b as u32 }

    /// Whether more bytes were read than the input has: the input was
    /// truncated (after decoding everything an [Encoder] wrote, the decoder
    /// has read exactly its bytes).
    pub fn overrun(&self) -> bool { self.pos > self.inp.len() }

    /// Decodes a bit coded with probability `p1` of a 1.
    pub fn bit(&mut self, p1: u32) -> u32 {
        let xmid = self.x1 + (((self.x2 - self.x1) as u64 * p1 as u64) >> 16) as u32;
        let bit = (self.x <= xmid) as u32;
        if bit == 1 { self.x2 = xmid } else { self.x1 = xmid + 1 }
        while (self.x1 ^ self.x2) & 0xFF00_0000 == 0 { self.x1 <<= 8; self.x2 = (self.x2 << 8) | 255; self.x = (self.x << 8) | self.byte(); }
        bit
    }
}

/// Either side of the coder: codes `bit` (ignored when decoding) with
/// probability `p1` (16-bit) of a 1 and returns the bit.
pub trait Coder {
    fn code(&mut self, bit: u32, p1: u32) -> u32;
    /// Whether a decoder has read past its input (never for an encoder), see
    /// [Decoder::overrun]; loops over stored counts stop on it.
    fn overrun(&self) -> bool { false }
}
impl Coder for Encoder { fn code(&mut self, bit: u32, p1: u32) -> u32 { self.bit(bit, p1); bit } }
impl Coder for Decoder<'_> {
    fn code(&mut self, _bit: u32, p1: u32) -> u32 { self.bit(p1) }
    fn overrun(&self) -> bool { Decoder::overrun(self) }
}

/// Codes a bit with an adaptive model and updates it.
#[inline]
pub fn code_bit<C: Coder>(
    c: &mut C,
    counter: &mut u32,
    bit: u32,
) -> u32 {
    let b = c.code(bit, p1(*counter));
    update(counter, b);
    b
}

/// Codes `v` (< 2^nbits) through an adaptive binary tree of `2^nbits`
/// models (`tree[1..]`, most significant bit first).
pub fn code_tree<C: Coder>(
    c: &mut C,
    tree: &mut [u32],
    nbits: u32,
    v: u32,
) -> u32 {
    let mut node = 1usize;
    let mut got = 0u32;
    for b in (0..nbits).rev() {
        let bit = code_bit(c, &mut tree[node], (v >> b) & 1);
        node = node * 2 + bit as usize;
        got = (got << 1) | bit;
    }
    got
}

/// Codes `v` as an Elias-gamma bucket (the bit length of v + 1, through a
/// 6-bit adaptive tree of 64 models) and the bits below the top one raw.
pub fn code_gamma<C: Coder>(
    c: &mut C,
    tree: &mut [u32],
    v: u64,
) -> u64 {
    let bucket = code_tree(c, tree, 6, 64 - (v + 1).leading_zeros() - 1) + 1;
    let mut x = 1u64;
    for b in (0..bucket - 1).rev() { x = (x << 1) | c.code((((v + 1) >> b) & 1) as u32, 32768) as u64; }
    x - 1
}

/// Adaptive Elias-gamma code of a number w ≥ 1: the bit length of w through a
/// 6-bit adaptive tree and the bits below its top one, the first
/// [GammaModel::MODELLED_BITS] adaptive given the bit length and the bits
/// above them, the rest raw.
pub struct GammaModel {
    bits: Vec<u32>,
    low: Vec<u32>,
}

impl Default for GammaModel {
    fn default() -> Self { Self::new() }
}

impl GammaModel {
    pub const MODELLED_BITS: u32 = 4;

    pub fn new() -> Self {
        GammaModel { bits: vec![C0; 64], low: vec![C0; 65 << Self::MODELLED_BITS] }
    }

    /// Codes (or decodes, `w` ignored) a number w ≥ 1.
    pub fn code<C: Coder>(
        &mut self,
        c: &mut C,
        w: u64,
    ) -> u64 {
        let n = code_tree(c, &mut self.bits, 6, 63 - w.max(1).leading_zeros()) + 1;
        let mut x = 1_u64;
        for (j, b) in (0..n - 1).rev().enumerate() {
            let truth = ((w >> b) & 1) as u32;
            let bit = if (j as u32) < Self::MODELLED_BITS {
                code_bit(c, &mut self.low[((n as usize) << Self::MODELLED_BITS) + x as usize], truth)
            } else {
                c.code(truth, 32768)
            };
            x = (x << 1) | bit as u64;
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rates_match_the_float_formula() {
        // As first written (and as the tig coder's files were made with).
        for n in 0..1024 { assert_eq!(DIV[n], (65536.0 / (n as f64 + 1.5)) as u32); }
    }

    #[test]
    fn roundtrip() {
        let values: Vec<u64> = vec![0, 1, 2, 3, 7, 8, 1000, 123_456_789, 0, 0, 5];
        let bits: Vec<u32> = (0..500).map(|i| (i % 7 == 0) as u32).collect();
        let mut enc = Encoder::new();
        let (mut tree, mut ctr) = (vec![C0; 64], C0);
        for &v in &values { code_gamma(&mut enc, &mut tree, v); }
        for &b in &bits { code_bit(&mut enc, &mut ctr, b); }
        let bytes = enc.finish();
        let mut dec = Decoder::new(&bytes);
        let (mut tree, mut ctr) = (vec![C0; 64], C0);
        let got: Vec<u64> = values.iter().map(|_| code_gamma(&mut dec, &mut tree, 0)).collect();
        let got_bits: Vec<u32> = bits.iter().map(|_| code_bit(&mut dec, &mut ctr, 0)).collect();
        assert_eq!(got, values);
        assert_eq!(got_bits, bits);
        assert!(!dec.overrun());
        // A truncated stream is noticed.
        let mut short = Decoder::new(&bytes[..bytes.len() - 5]);
        let (mut tree, mut ctr) = (vec![C0; 64], C0);
        for _ in &values { code_gamma(&mut short, &mut tree, 0); }
        for _ in &bits { code_bit(&mut short, &mut ctr, 0); }
        assert!(short.overrun());
        // 500 bits with a 1-in-7 skew cost ~0.6 bits each (~37 bytes) plus the
        // values (~13 bytes), well below 500 / 8 + 13.
        assert!(bytes.len() < 60, "{} bytes", bytes.len());
    }
}
