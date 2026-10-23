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
//! Logistic mixing: the stretch (ln(p / (1 - p))) and squash (its inverse)
//! tables shared by the context-mixing models, and the adaptive probability
//! map ([Apm]) that refines their mixed predictions.
//!
//! The tables are computed in fixed-point integer arithmetic (no
//! floating-point logarithms or exponentials), so they are the same on every
//! platform and the models that use them predict exactly the same
//! probabilities in the encoder and the decoder.

use super::arith::update;

/// Stretch range: +-SMAX in 1/256 nat.
pub const SMAX: i32 = 4095;

// Fixed-point arithmetic for the tables (Q62 in u128 / i128) ----------------------

const Q: u32 = 62;
const QONE: u128 = 1 << Q;
/// ln 2 in Q62, rounded.
const LN2_Q: u128 = 3_196_577_161_300_663_915;

/// Q62 product (both factors below 2^65).
fn q_mul(a: u128, b: u128) -> u128 { (a * b) >> Q }

/// exp(-1/256) in Q62, by its Taylor series.
fn q_exp_m1_256() -> u128 {
    let (mut sum, mut term): (i128, i128) = (QONE as i128, QONE as i128);
    for n in 1.. {
        term = -term / (256 * n);
        if term == 0 { break }
        sum += term;
    }
    sum as u128
}

/// exp(-x / 256) in Q62 for x >= 0, by repeated squaring of exp(-1/256).
fn q_exp_neg(x: u64, e1: u128) -> u128 {
    let (mut result, mut base, mut e) = (QONE, e1, x);
    while e > 0 {
        if e & 1 == 1 { result = q_mul(result, base); }
        base = q_mul(base, base);
        e >>= 1;
    }
    result
}

/// log2(n) in Q62 for n >= 1: the integer part from the bit length, the
/// fraction bit by bit by squaring the mantissa.
fn q_log2(n: u64) -> u128 {
    let k = 63 - n.leading_zeros();
    let mut m: u128 = ((n as u128) << Q) >> k; // in [1, 2)
    let mut frac: u128 = 0;
    for b in (0..Q).rev() {
        m = q_mul(m, m);
        if m >= 2 * QONE { m >>= 1; frac |= 1 << b; }
    }
    ((k as u128) << Q) | frac
}

/// A signed Q62 value rounded to the nearest integer, half away from zero.
fn q_round(v: i128) -> i64 {
    let half = (QONE / 2) as i128;
    if v >= 0 { ((v + half) >> Q) as i64 } else { -(((-v + half) >> Q) as i64) }
}

/// The stretch and squash tables, computed in integer arithmetic so that every
/// platform gets the same values (each is the correctly rounded value: none
/// lies within 3.9e-6 of a rounding tie, far above the Q62 error).
pub struct Tables { pub stretch: Vec<i32>, pub squash: Vec<u32> }
impl Tables {
    pub fn new() -> Tables {
        // stretch[i] = round(256 ln(p / (1 - p))) within +-SMAX, for the 16-bit
        // probability p = (i + 0.5) / 65536 = (2i + 1) / 131072.
        let (ln2_hi, ln2_lo) = ((LN2_Q >> 31) as i128, (LN2_Q & ((1 << 31) - 1)) as i128);
        let stretch = (0..65536u64).map(|i| {
            let log2 = q_log2(2 * i + 1) as i128 - q_log2(131_071 - 2 * i) as i128;
            // 256 log2 ln 2, with ln 2 split in two halves to stay within 128 bits.
            let v = log2 * 256;
            let ln = ((v * ln2_hi) >> (Q - 31)) + ((v * ln2_lo) >> Q);
            q_round(ln).clamp(-SMAX as i64, SMAX as i64) as i32
        }).collect();
        // squash[j] = round(65536 / (1 + exp(-x / 256))) within 1..=65535, for
        // x = j - SMAX.
        let e1 = q_exp_m1_256();
        let squash = (0..2 * SMAX + 1).map(|j| {
            let x = j - SMAX;
            let e = q_exp_neg(x.unsigned_abs() as u64, e1); // exp(-|x| / 256)
            let (num, den) = if x >= 0 { (65536 << Q, QONE + e) } else { (65536 * e, QONE + e) };
            ((num + den / 2) / den).clamp(1, 65535) as u32
        }).collect();
        Tables { stretch, squash }
    }
    /// Stretch of a counter (probability in its high 22 bits).
    #[inline] pub fn st(&self, c: u32) -> i32 { self.stretch[(c >> 16) as usize] }
    /// 16-bit probability of a stretch.
    #[inline] pub fn sq(&self, x: i32) -> u32 { self.squash[(x.clamp(-SMAX, SMAX) + SMAX) as usize] }
}

/// Bins of an [Apm] per context: the clamped stretch -2047..=2047 in steps
/// of 128.
const APM_BINS: usize = 33;

/// How an [Apm] bin moves towards a coded bit.
#[derive(Clone, Copy, Debug)]
pub enum ApmRate {
    /// As an adaptive bit model ([update]): fast at first, slower as its
    /// count grows.
    Adaptive,
    /// By a fixed fraction (in 1/65536) of the distance to the bit.
    Fixed(i64),
}

/// An adaptive probability map: refines a prediction (a stretch) by a
/// context. Each context has [APM_BINS] bins over the stretch, starting at the
/// identity (the squash of the bin's stretch); a prediction interpolates
/// between the two bins around it, and only the nearer one learns the coded
/// bit. The bins are counters of adaptive bit models (probability in the high
/// 22 bits); a [ApmRate::Fixed] map keeps their counts at 0.
pub struct Apm {
    bins: Vec<u32>,
    rate: ApmRate,
}

impl Apm {
    /// A map with `n_contexts` contexts.
    pub fn new(
        t: &Tables,
        n_contexts: usize,
        rate: ApmRate,
    ) -> Apm {
        let bins = (0..n_contexts * APM_BINS).map(|j| {
            let x = ((j % APM_BINS) as i32 - 16) * 128;
            (t.sq(x) << 6).min((1 << 22) - 1) << 10
        }).collect();
        Apm { bins, rate }
    }

    /// The refined 16-bit probability of stretch `s` in context `ctx`, and the
    /// bin that learns the coded bit ([Apm::update]).
    #[inline]
    pub fn refine(
        &self,
        ctx: usize,
        s: i32,
    ) -> (u32, usize) {
        let sc = s.clamp(-2047, 2047) + 2048;
        let base = ctx * APM_BINS + (sc >> 7) as usize;
        let fr = sc & 127;
        let p = ((self.bins[base] >> 16) as i64 * (128 - fr) as i64 + (self.bins[base + 1] >> 16) as i64 * fr as i64) >> 7;
        (p as u32, if fr < 64 { base } else { base + 1 })
    }

    /// Moves `bin` towards the coded `bit`.
    #[inline]
    pub fn update(
        &mut self,
        bin: usize,
        bit: u32,
    ) {
        let counter = &mut self.bins[bin];
        match self.rate {
            ApmRate::Adaptive => update(counter, bit),
            ApmRate::Fixed(rate) => {
                let p = (*counter >> 10) as i64;
                let p = p + ((((bit as i64) << 22) - p) * rate >> 16);
                *counter = (p as u32) << 10;
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fnv1a(bytes: impl Iterator<Item = u8>) -> u64 {
        bytes.fold(0xcbf2_9ce4_8422_2325, |h, b| (h ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3))
    }

    #[test]
    fn tables_are_fixed() {
        // The tables define the format: checksums of their little-endian i16 /
        // u16 values (as first written with f64 ln and exp).
        let t = Tables::new();
        assert_eq!(fnv1a(t.stretch.iter().flat_map(|&x| (x as i16).to_le_bytes())), 0xbccd_3940_e4c0_fead);
        assert_eq!(fnv1a(t.squash.iter().flat_map(|&x| (x as u16).to_le_bytes())), 0x0142_3c05_9f88_cf45);
        assert_eq!((t.stretch[0], t.stretch[32767], t.stretch[32768], t.stretch[65535]), (-3017, 0, 0, 3017));
        assert_eq!((t.squash[0], t.squash[4095], t.squash[8190]), (1, 32768, 65535));
    }

    #[test]
    fn tables_are_monotonic() {
        let t = Tables::new();
        assert_eq!(t.stretch.len(), 65536);
        assert_eq!(t.squash.len(), 2 * SMAX as usize + 1);
        assert!(t.stretch.windows(2).all(|w| w[0] <= w[1]));
        assert!(t.squash.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!(t.sq(0), 32768);
    }

    #[test]
    fn apm_matches_its_first_versions() {
        // The choice model's APM as first written: 22-bit probabilities in
        // i32 with a fixed rate; the tig coder's: adaptive counters.
        let t = Tables::new();
        let (n_contexts, rate) = (3, 1311);
        let mut fixed = Apm::new(&t, n_contexts, ApmRate::Fixed(rate));
        let mut adaptive = Apm::new(&t, n_contexts, ApmRate::Adaptive);
        let mut old_fixed: Vec<i32> = (0..n_contexts * 33).map(|i| (t.sq(((i % 33) as i32 - 16) * 128) as i32) << 6).collect();
        let mut old_adaptive: Vec<u32> = (0..n_contexts * 33).map(|j| (t.sq(((j % 33) as i32 - 16) * 128) << 6).min((1 << 22) - 1) << 10).collect();
        let mut x = 12345_u64;
        for _ in 0..100_000 {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            let (ctx, s, bit) = ((x % 3) as usize, ((x >> 8) % 8191) as i32 - 4095, ((x >> 32) % 10 < 7) as u32);
            let sc = s.clamp(-2047, 2047) + 2048;
            let (base, fr) = (ctx * 33 + (sc >> 7) as usize, sc & 127);
            let near = if fr < 64 { base } else { base + 1 };

            let p = ((old_fixed[base] >> 6) as i64 * (128 - fr) as i64 + (old_fixed[base + 1] >> 6) as i64 * fr as i64) >> 7;
            assert_eq!(fixed.refine(ctx, s), (p as u32, near));
            fixed.update(near, bit);
            old_fixed[near] += ((((bit as i64) << 22) - old_fixed[near] as i64) * rate >> 16) as i32;

            let p = ((old_adaptive[base] >> 16) as i64 * (128 - fr) as i64 + (old_adaptive[base + 1] >> 16) as i64 * fr as i64) >> 7;
            assert_eq!(adaptive.refine(ctx, s), (p as u32, near));
            adaptive.update(near, bit);
            update(&mut old_adaptive[near], bit);
        }
    }
}
