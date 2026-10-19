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
//! DNA context-mixing coder for the index tigs.
//!
//! The tigs are a spectrum-preserving string set of a large, redundant
//! _k_-mer set: most bases repeat earlier ones with substitutions, and every
//! _k_-mer occurs once. [encode] codes them with an adaptive model in integer
//! arithmetic (so that [decode] repeats its predictions exactly) driving a
//! binary arithmetic coder (INDEX_JOURNAL.md, "Index tigs: a DNA-specific
//! coder").
//!
//! **Layout.** The long tigs (at least 2_k_ bases) come first, longest first.
//! Then rounds of anchored short tigs: a short tig whose first or last _k_ - 1
//! bases occur in the previous region (round 1: the long tigs; round _r_: the
//! tigs of round _r_ - 1), on either strand, is written with that end first,
//! and each round is sorted by anchor position. Anchors on the forward strand
//! of the region are preferred. Then the short tigs with no anchor, longest
//! first. [decode] returns the tigs in this order and orientation: the same
//! _k_-mers, so the same SBWT.
//!
//! **Format.** Magic `TAC5`, _k_ (u16), options (u8: bit 0 end flags, always
//! set; bit 1 reverse-strand anchors, set by [encode]), number of
//! rounds (u8), total bases (u64), numbers of long tigs, of unanchored short
//! tigs and of tigs in each round (u32 each, little-endian), then one
//! arithmetic-coded stream:
//! - the lengths of the long and of the unanchored tigs: the first raw, then
//!   the differences from the previous one (adaptive Elias-gamma bucket plus
//!   raw low bits);
//! - per round, the anchor positions as differences, likewise, each followed
//!   by its strand (an adaptive bit) if reverse-strand anchors are on;
//! - the bases. An anchored tig's first _k_ - 1 bases are copied from its
//!   anchor, reverse complemented for a reverse-strand anchor (the model
//!   still sees them); from base _k_ on, an end flag precedes each of its
//!   bases.
//!
//! **Model**, per base as two binary decisions (high bit, then low bit):
//! - order-_n_ context counters for _n_ in [ORDERS] (direct tables while they
//!   fit, hashed with a check above);
//! - copy models: [N_FWD] forward and [N_RC] reverse-complement ones, each
//!   following a different one of the latest positions of the hashed
//!   [MINLEN]-mer, plus a [SHORT]-mer one looked up only in the first
//!   [SHORT_GATE] bases of a tig; tolerant to substitutions; after a miss a
//!   longer candidate replaces the match; all are dropped at a tig start;
//! - _k_-mer exclusion: after an exact _k_ - 1 base match within one tig the
//!   base that followed there would repeat a _k_-mer, so it is impossible;
//! - two-layer mixing in the logistic domain: [N_MIX] mixers whose weight
//!   sets are selected per node by the copy state (longest match, agreement
//!   of the copy models, tig start), the order-2 context, the order-16 count,
//!   and the base the longest copy predicts (with the agreement, or with a
//!   differing predicted base); a final mixer over their outputs;
//! - two APMs: by copy state, and by predicted base and agreement.
//!
//! The stretch and squash tables come from [super::mix], computed in integer
//! arithmetic, so they are the same on every platform.
use std::collections::HashMap;
use std::io::{Read, Write};

use super::arith::{code_gamma, update, Coder, Decoder, Encoder, C0};
use super::mix::{Tables, SMAX};

type E = Box<dyn std::error::Error>;

const MAGIC: &[u8; 4] = b"TAC5";

/// Orders of the context models.
const ORDERS: [usize; 13] = [1, 2, 3, 4, 6, 8, 11, 12, 14, 16, 18, 20, 24];
/// Index of order 16 in [ORDERS].
const ORDER16: usize = 9;
/// Count-weighted inputs: the highest orders, from this index of [ORDERS] on.
const CONF_FROM: usize = 8;
const N_CONF: usize = ORDERS.len() - CONF_FROM;
/// Largest table sizes (log2 of the number of contexts or buckets).
const MAX_HASH_BITS: u32 = 25;
const MAX_MATCH_BITS: u32 = 26;
/// Copy models: minimum match, forward and reverse-complement candidates.
const MINLEN: usize = 20;
const N_FWD: usize = 4;
const N_RC: usize = 4;
/// The short copy model, looked up only in the first SHORT_GATE bases of a tig.
const SHORT: usize = 14;
const SHORT_GATE: usize = 40;
const N_MODELS: usize = N_FWD + N_RC + 1;
/// Weight sets of the first mixer and of the first APM, per node.
const N_SEL: usize = 7 * 3 * 2;
const N_IN: usize = ORDERS.len() + N_MODELS + N_CONF + 1;
/// First-layer mixers, by weight-set selector: copy state; order-2 context;
/// order-16 count x agreement; predicted base x agreement x longest match;
/// predicted base x a differing predicted base x longest match.
const N_MIX: usize = 5;
const MIX_SETS: [usize; N_MIX] = [N_SEL, 16, 7 * 3, 5 * 3 * 7, 5 * 5 * 7];
/// Mixer weight 1.0 and learning rate (dw = x * err >> LR_SHIFT).
const W_ONE: i64 = 1 << 20;
const LR_SHIFT: u32 = 13;
/// Anchoring rounds at most.
const MAX_ROUNDS: usize = 4;

fn hash64(x: u64, n: usize) -> u64 { (x ^ (n as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)).wrapping_mul(0xD6E8_FEB8_6659_FD93).rotate_left(29).wrapping_mul(0x9E37_79B9_7F4A_7C15) }

#[derive(Clone, Copy)]
struct Match { ptr: usize, exact: usize, hist: u32, active: bool }

/// A copy model: follows candidate `slot` of hash table `table` on one strand.
#[derive(Clone, Copy)]
struct CopyModel { rc: bool, table: usize, slot: usize, m: Match }
impl CopyModel {
    fn predicted(&self, text: &[u8]) -> Option<u8> {
        if !self.m.active { return None }
        Some(if self.rc { 3 - text[self.m.ptr] } else { text[self.m.ptr] })
    }
    /// Match length bucket; 6 for a full k - 1 base match.
    fn len_bucket(&self, k: usize) -> usize {
        if !self.m.active { return 0 }
        if self.m.exact >= k - 1 { return 6 }
        match self.m.exact { 0 => 1, 1..=15 => 2, 16..=31 => 3, 32..=63 => 4, _ => 5 }
    }
    fn bucket(&self, k: usize) -> usize { self.len_bucket(k) + 7 * (self.m.hist & 0xFFFF).count_ones().min(3) as usize }
}

/// Positions after the latest occurrences of each hashed `minlen`-mer,
/// `slots` per bucket, newest first.
struct HTab { minlen: usize, slots: usize, bits: u32, data: Vec<u32> }
impl HTab {
    fn bucket(&self, key: u64) -> usize { (hash64(key, 99 + self.minlen) >> (64 - self.bits)) as usize * self.slots }
}

struct Model {
    k: usize,
    t: Tables,
    hash_bits: u32,
    n_direct: usize,
    direct: Vec<Vec<u32>>,
    hashed: Vec<Vec<u32>>,
    tabs: Vec<HTab>,
    models: Vec<CopyModel>,
    mprob: Vec<u32>,
    weights: Vec<Vec<i32>>,
    fweights: Vec<i32>,
    apm: Vec<u32>,
    apm2: Vec<u32>,
    f: u64,
    r: u64,
    // Per-base state between predict and update.
    slots: [usize; ORDERS.len()],
    inputs: [[i32; N_IN]; 3],
    msel: [[usize; N_MIX]; 3],
    psub: [[u32; N_MIX]; 3],
    finputs: [[i32; N_MIX + 1]; 3],
    pmix: [u32; 3],
    apm_at: [(usize, i32); 3],
    apm2_at: [usize; 3],
    sel: usize,
    preds: [Option<u8>; N_MODELS],
}

impl Model {
    /// A model for `k`-mers and about `n` bases (the table sizes follow `n`).
    fn new(k: usize, n: usize) -> Model {
        let lg = usize::BITS - n.max(1024).saturating_sub(1).leading_zeros(); // ceil(log2(n))
        let hash_bits = lg.saturating_sub(3).clamp(10, MAX_HASH_BITS);
        let match_bits = lg.saturating_sub(2).clamp(10, MAX_MATCH_BITS);
        let n_direct = ORDERS.iter().take_while(|&&o| 2 * o as u32 <= hash_bits).count();
        let t = Tables::new();
        let apm_init = |n: usize| -> Vec<u32> { (0..n * 33).map(|j| { let x = ((j % 33) as i32 - 16) * 128; (t.sq(x) << 6).min((1 << 22) - 1) << 10 }).collect() };
        let (apm, apm2) = (apm_init(3 * N_SEL), apm_init(3 * 16));
        let main_slots = N_FWD.max(N_RC);
        let tabs = vec![
            HTab { minlen: MINLEN, slots: main_slots, bits: match_bits, data: vec![u32::MAX; (1usize << match_bits) * main_slots] },
            HTab { minlen: SHORT, slots: 1, bits: match_bits, data: vec![u32::MAX; 1usize << match_bits] },
        ];
        let m0 = Match { ptr: 0, exact: 0, hist: 0, active: false };
        let mut models = Vec::new();
        for s in 0..N_FWD { models.push(CopyModel { rc: false, table: 0, slot: s, m: m0 }); }
        for s in 0..N_RC { models.push(CopyModel { rc: true, table: 0, slot: s, m: m0 }); }
        models.push(CopyModel { rc: false, table: 1, slot: 0, m: m0 });
        Model {
            k, t,
            hash_bits, n_direct,
            direct: ORDERS[..n_direct].iter().map(|&o| vec![C0; 4usize.pow(o as u32) * 4]).collect(),
            hashed: ORDERS[n_direct..].iter().map(|_| vec![0u32; (1usize << hash_bits) * 4]).collect(),
            tabs, models,
            mprob: vec![C0; N_MODELS * 28 * 2],
            weights: MIX_SETS.iter().map(|&n| vec![(W_ONE / 4) as i32; 3 * n * N_IN]).collect(),
            fweights: (0..3 * N_SEL * (N_MIX + 1)).map(|j| if j % (N_MIX + 1) == N_MIX { 0 } else { (W_ONE / N_MIX as i64) as i32 }).collect(),
            apm, apm2,
            f: 0, r: 0,
            slots: [0; ORDERS.len()], inputs: [[0; N_IN]; 3], msel: [[0; N_MIX]; 3], psub: [[0; N_MIX]; 3], finputs: [[0; N_MIX + 1]; 3], pmix: [0; 3], apm_at: [(0, 0); 3], apm2_at: [0; 3], sel: 0, preds: [None; N_MODELS],
        }
    }

    fn slot(&mut self, oi: usize) -> usize {
        let o = ORDERS[oi];
        let ctx = self.f & ((1u64 << (2 * o)) - 1);
        if oi < self.n_direct { return (ctx as usize) * 4 }
        let h = hash64(ctx, o);
        let i = ((h >> (64 - self.hash_bits)) as usize) * 4;
        let chk = (h & 0xFFFF) as u32 | 1 << 16;
        let t = &mut self.hashed[oi - self.n_direct];
        if t[i + 3] != chk { t[i] = C0; t[i + 1] = C0; t[i + 2] = C0; t[i + 3] = chk; }
        i
    }
    fn counter(&self, oi: usize, node: usize) -> u32 {
        if oi < self.n_direct { self.direct[oi][self.slots[oi] + node] } else { self.hashed[oi - self.n_direct][self.slots[oi] + node] }
    }
    fn longest_bucket(&self) -> usize { self.models.iter().map(|cm| cm.len_bucket(self.k)).max().unwrap() }

    /// Probabilities (16-bit) of a 1 at the three nodes for the base at `i`,
    /// and the excluded bases (bit mask). `text[..i]` is known.
    fn predict(&mut self, text: &[u8], i: usize, pos_in_tig: usize, tig_of: &dyn Fn(usize) -> usize) -> ([u32; 3], u8) {
        let k = self.k;
        if pos_in_tig == 0 { for cm in self.models.iter_mut() { cm.m.active = false; } }
        for oi in 0..ORDERS.len() { self.slots[oi] = self.slot(oi); }
        for j in 0..N_MODELS { self.preds[j] = self.models[j].predicted(text); }
        // Longest match bucket, agreement (none / agree / disagree), tig start.
        let lb = self.longest_bucket();
        let mut agree = 0;
        let mut first: Option<u8> = None;
        let mut n_pred = 0;
        for p in self.preds.iter().flatten() { n_pred += 1; match first { None => { first = Some(*p); agree = 1; }, Some(q) if q != *p => { agree = 2; }, _ => {} } }
        if agree == 1 && n_pred < 2 { agree = 0; }
        self.sel = (lb * 3 + agree) * 2 + (pos_in_tig < 40) as usize;
        let ctx2 = (self.f & 15) as usize;
        // The longest copy's predicted base (+1, 0: none) and a differing one.
        let lp = (0..N_MODELS).filter(|&j| self.models[j].m.active).max_by_key(|&j| self.models[j].m.exact).and_then(|j| self.preds[j]).map_or(0, |e| e as usize + 1);
        let alt = if lp == 0 { 0 } else { self.preds.iter().flatten().find(|&&e| e as usize + 1 != lp).map_or(0, |&e| e as usize + 1) };
        let key2 = lp * 3 + agree;
        let mut p = [0u32; 3];
        for node in 0..3 {
            let mut x = [0i32; N_IN];
            for oi in 0..ORDERS.len() { x[oi] = self.t.st(self.counter(oi, node)); }
            for j in 0..N_MODELS {
                x[ORDERS.len() + j] = match self.preds[j] {
                    Some(e) if node == 0 || (node - 1) as u8 == e >> 1 => {
                        let eb = if node == 0 { e >> 1 } else { e & 1 };
                        let s = self.t.st(self.mprob[(j * 28 + self.models[j].bucket(k)) * 2 + (node > 0) as usize]);
                        if eb == 1 { s } else { -s }
                    },
                    _ => 0,
                };
            }
            for q in 0..N_CONF {
                let c = self.counter(CONF_FROM + q, node);
                let n = (c & 1023) as i32;
                x[ORDERS.len() + N_MODELS + q] = self.t.st(c) * n / (n + 2);
            }
            x[N_IN - 1] = 77;
            let n16 = self.counter(ORDER16, node) & 1023;
            let cb = match n16 { 0 => 0, 1 => 1, 2 => 2, 3..=4 => 3, 5..=8 => 4, 9..=16 => 5, _ => 6 };
            let sels = [self.sel, ctx2, cb * 3 + agree, (lp * 3 + agree) * 7 + lb, (lp * 5 + alt) * 7 + lb];
            let mut fx = [0i32; N_MIX + 1];
            for mj in 0..N_MIX {
                let w = &self.weights[mj][(node * MIX_SETS[mj] + sels[mj]) * N_IN..][..N_IN];
                let dot: i64 = x.iter().zip(w).map(|(&a, &b)| a as i64 * b as i64).sum();
                let sj = ((dot / W_ONE) as i32).clamp(-SMAX, SMAX);
                self.psub[node][mj] = self.t.sq(sj);
                fx[mj] = sj;
            }
            fx[N_MIX] = 77;
            self.msel[node] = sels;
            let w = &self.fweights[(node * N_SEL + self.sel) * (N_MIX + 1)..][..N_MIX + 1];
            let dot: i64 = fx.iter().zip(w).map(|(&a, &b)| a as i64 * b as i64).sum();
            let s = (dot / W_ONE) as i32;
            let pm = self.t.sq(s);
            self.finputs[node] = fx;
            self.pmix[node] = pm;
            self.inputs[node] = x;
            // APMs over 33 bins of the clamped stretch.
            let sc = s.clamp(-2047, 2047) + 2048;
            let j = (sc >> 7) as usize; let fr = sc & 127;
            let base = (node * N_SEL + self.sel) * 33 + j;
            self.apm_at[node] = (base, fr);
            let pa = ((self.apm[base] >> 16) as i64 * (128 - fr) as i64 + (self.apm[base + 1] >> 16) as i64 * fr as i64) >> 7;
            let b2 = (node * 16 + key2) * 33 + j;
            self.apm2_at[node] = b2;
            let pa2 = ((self.apm2[b2] >> 16) as i64 * (128 - fr) as i64 + (self.apm2[b2 + 1] >> 16) as i64 * fr as i64) >> 7;
            p[node] = ((pm as i64 + pa + 2 * pa2) / 4).clamp(1, 65535) as u32;
        }
        // k-mer exclusion.
        let mut excluded = 0u8;
        if pos_in_tig >= k - 1 {
            for j in 0..N_MODELS {
                let cm = self.models[j];
                if let (Some(e), true) = (self.preds[j], cm.m.exact >= k - 1) {
                    let q = cm.m.ptr;
                    let ok = if cm.rc { q + k - 1 < i && tig_of(q) == tig_of(q + k - 1) } else { q >= k - 1 && tig_of(q) == tig_of(q - (k - 1)) };
                    if ok { excluded |= 1 << e; }
                }
            }
        }
        (p, excluded)
    }

    fn update(&mut self, text: &[u8], i: usize, pos_in_tig: usize, sym: u8) {
        let k = self.k;
        let b1 = (sym >> 1) as u32;
        let b0 = (sym & 1) as u32;
        let round = 1i64 << (LR_SHIFT - 1);
        for (node, bit) in [(0usize, b1), (1 + b1 as usize, b0)] {
            for mj in 0..N_MIX {
                let err = ((bit as i64) << 16) - self.psub[node][mj] as i64;
                let w = &mut self.weights[mj][(node * MIX_SETS[mj] + self.msel[node][mj]) * N_IN..][..N_IN];
                for (wi, &x) in w.iter_mut().zip(&self.inputs[node]) { *wi += ((x as i64 * err + round) >> LR_SHIFT) as i32; }
            }
            let err = ((bit as i64) << 16) - self.pmix[node] as i64;
            let w = &mut self.fweights[(node * N_SEL + self.sel) * (N_MIX + 1)..][..N_MIX + 1];
            for (wi, &x) in w.iter_mut().zip(&self.finputs[node]) { *wi += ((x as i64 * err + round) >> LR_SHIFT) as i32; }
            let (b, fr) = self.apm_at[node];
            let (a, a2) = if fr < 64 { (b, self.apm2_at[node]) } else { (b + 1, self.apm2_at[node] + 1) };
            update(&mut self.apm[a], bit);
            update(&mut self.apm2[a2], bit);
            for oi in 0..ORDERS.len() {
                let s = self.slots[oi] + node;
                if oi < self.n_direct { update(&mut self.direct[oi][s], bit) } else { update(&mut self.hashed[oi - self.n_direct][s], bit) }
            }
            for j in 0..N_MODELS {
                if let Some(e) = self.preds[j] {
                    if node == 0 || (node - 1) as u8 == e >> 1 {
                        let eb = if node == 0 { e >> 1 } else { e & 1 } as u32;
                        update(&mut self.mprob[(j * 28 + self.models[j].bucket(k)) * 2 + (node > 0) as usize], (eb == bit) as u32);
                    }
                }
            }
        }
        self.f = (self.f << 2) | sym as u64;
        self.r = (self.r >> 2) | ((3 - sym as u64) << 62);
        for j in 0..N_MODELS {
            if let Some(e) = self.preds[j] {
                let cm = &mut self.models[j];
                let m = &mut cm.m;
                let hit = e == sym;
                m.hist = (m.hist << 1) | (!hit) as u32;
                if hit { m.exact += 1; } else { m.exact = 0; }
                let next = if cm.rc { m.ptr.checked_sub(1) } else { Some(m.ptr + 1).filter(|&p| p <= i) };
                let drop = !hit && (m.hist & 0xFFFF).count_ones() > 8;
                match next { Some(p) if !drop => m.ptr = p, _ => { m.active = false; } }
            }
        }
        // Lookups: each copy model takes its slot's candidate if it beats its current match.
        let mask = |o: usize| (1u64 << (2 * o)) - 1;
        let buckets: Vec<Option<(usize, usize)>> = self.tabs.iter().map(|t| if i + 1 >= t.minlen { Some((t.bucket(self.f & mask(t.minlen)), t.bucket(self.r >> (2 * (32 - t.minlen))))) } else { None }).collect();
        for j in 0..N_MODELS {
            let (rc, ti, slot) = (self.models[j].rc, self.models[j].table, self.models[j].slot);
            let Some((bf, br)) = buckets[ti] else { continue };
            if ti == 1 && pos_in_tig + 1 >= SHORT_GATE { continue }
            let m = self.models[j].m;
            if m.active && m.hist & 1 == 0 { continue }
            let t = &self.tabs[ti];
            let minlen = t.minlen;
            let cand = t.data[(if rc { br } else { bf }) + slot];
            if cand == u32::MAX { continue }
            let q = cand as usize;
            let (ptr, l) = if !rc {
                if q > i { continue }
                let mut l = 0;
                while l < 400 && l < q && text[q - 1 - l] == text[i - l] { l += 1; }
                (q, l)
            } else {
                if q <= minlen { continue }
                let start = q - minlen;
                let mut l = 0;
                while l < 400 && start + l < i + 1 - l && text[start + l] == 3 - text[i - l] { l += 1; }
                (start - 1, l)
            };
            if l < minlen || (m.active && l <= m.exact) { continue }
            if self.models.iter().enumerate().any(|(o, cm)| o != j && cm.rc == rc && cm.m.active && cm.m.ptr == ptr) { continue }
            self.models[j].m = Match { ptr, exact: l, hist: 0, active: true };
        }
        for (ti, t) in self.tabs.iter_mut().enumerate() {
            if let Some((bf, _)) = buckets[ti] {
                let d = &mut t.data[bf..bf + t.slots];
                d.copy_within(0..t.slots - 1, 1);
                d[0] = (i + 1) as u32;
            }
        }
    }
}

/// Probabilities of the four bases with the excluded ones (bit mask) removed, and their sum.
fn node_probs(p: [u32; 3], ex: u8) -> ([u64; 4], u64) {
    let mut q = [
        (65536 - p[0] as u64) * (65536 - p[1] as u64),
        (65536 - p[0] as u64) * p[1] as u64,
        p[0] as u64 * (65536 - p[2] as u64),
        p[0] as u64 * p[2] as u64,
    ];
    for (e, x) in q.iter_mut().enumerate() { if ex >> e & 1 == 1 { *x = 0; } }
    (q, q.iter().sum())
}
fn prob16(num: u64, den: u64) -> u32 { (((num as u128) << 16) / den as u128).clamp(1, 65535) as u32 }

/// Codes (or decodes) one base given the node probabilities and the excluded bases.
fn code_base<C: Coder>(c: &mut C, p: [u32; 3], ex: u8, sym: u8) -> u8 {
    if ex == 0 {
        let b1 = c.code((sym >> 1) as u32, p[0]);
        let b0 = c.code((sym & 1) as u32, p[1 + b1 as usize]);
        return ((b1 << 1) | b0) as u8
    }
    let (q, tot) = node_probs(p, ex);
    let hi = q[2] + q[3];
    // A bit is not coded if only one value remains.
    let b1 = if hi == 0 { 0 } else if hi == tot { 1 } else { c.code((sym >> 1) as u32, prob16(hi, tot)) };
    let (a, b) = (q[2 * b1 as usize], q[2 * b1 as usize + 1]);
    let b0 = if b == 0 { 0 } else if a == 0 { 1 } else { c.code((sym & 1) as u32, prob16(b, a + b)) };
    ((b1 << 1) | b0) as u8
}

/// Codes `n` lengths, longest first: the first as 32 raw bits, then the
/// differences from the previous one. Decoding appends them to `lens`.
fn code_lengths<C: Coder>(c: &mut C, lens: &mut Vec<usize>, n: usize) {
    let mut tree = vec![C0; 128];
    let mut prev = 0usize;
    for j in 0..n {
        let v = if j == 0 {
            let mut v = 0u64;
            for b in (0..32).rev() { v |= (c.code((lens.first().copied().unwrap_or(0) as u64 >> b & 1) as u32, 32768) as u64) << b; }
            v as usize
        } else {
            let d = if lens.len() > j { (prev - lens[j]) as u64 } else { 0 };
            prev - code_gamma(c, &mut tree, d) as usize
        };
        if lens.len() <= j { lens.push(v); }
        prev = v;
    }
}

// Layout -------------------------------------------------------------------------

const B1: u64 = 0x9E37_79B9_7F4A_7C15 | 1;
const B2: u64 = 0xC2B2_AE3D_27D4_EB4F | 1;
/// Hash of a sequence of base codes.
fn khash(s: &[u8]) -> u128 {
    let (mut a, mut b) = (0u64, 0u64);
    for &c in s { a = a.wrapping_mul(B1).wrapping_add(c as u64 + 1); b = b.wrapping_mul(B2).wrapping_add(c as u64 + 1); }
    ((a as u128) << 64) | b as u128
}
fn rc_codes(s: &[u8]) -> Vec<u8> { s.iter().rev().map(|&c| 3 - c).collect() }

/// Adds to `want` the positions (`region_start` + offset) of its `w`-mers in
/// the tigs of a region.
fn find_kmers(region: &[&[u8]], region_start: usize, w: usize, want: &mut HashMap<u128, Vec<usize>>) {
    let (p1, p2) = ((0..w - 1).fold(1u64, |x, _| x.wrapping_mul(B1)), (0..w - 1).fold(1u64, |x, _| x.wrapping_mul(B2)));
    let mut start = region_start;
    for t in region {
        let (mut a, mut b) = (0u64, 0u64);
        for (i, &c) in t.iter().enumerate() {
            if i >= w { let o = t[i - w] as u64 + 1; a = a.wrapping_sub(o.wrapping_mul(p1)); b = b.wrapping_sub(o.wrapping_mul(p2)); }
            a = a.wrapping_mul(B1).wrapping_add(c as u64 + 1); b = b.wrapping_mul(B2).wrapping_add(c as u64 + 1);
            if i + 1 >= w { if let Some(v) = want.get_mut(&(((a as u128) << 64) | b as u128)) { v.push(start + i + 1 - w); } }
        }
        start += t.len();
    }
}

/// Options byte of the format: the anchored tigs' ends are coded as flags
/// (always set); anchors may be on the reverse strand.
const OPT_END_FLAGS: u8 = 1;
const OPT_RC_ANCHORS: u8 = 2;

/// The tigs in coding order, see the module docs; an anchored tig is
/// (anchor position, on the reverse strand, sequence).
struct Plan { long: Vec<Vec<u8>>, rounds: Vec<Vec<(usize, bool, Vec<u8>)>>, rest: Vec<Vec<u8>> }

fn plan(tigs: Vec<Vec<u8>>, k: usize) -> Plan {
    let w = k - 1;
    let (long, mut rest): (Vec<Vec<u8>>, Vec<Vec<u8>>) = tigs.into_iter().partition(|t| t.len() >= 2 * k);
    let mut rounds: Vec<Vec<(usize, bool, Vec<u8>)>> = Vec::new();
    let mut region_start = 0usize;
    let mut region_len: usize = long.iter().map(|t| t.len()).sum();
    while rounds.len() < MAX_ROUNDS {
        let region: Vec<&[u8]> = match rounds.last() { None => long.iter().map(|t| t.as_slice()).collect(), Some(r) => r.iter().map(|(_, _, t)| t.as_slice()).collect() };
        let mut want: HashMap<u128, Vec<usize>> = HashMap::new();
        // Both ends of each tig, both strands.
        for t in &rest {
            let (head, tail) = (&t[..w], &t[t.len() - w..]);
            for key in [khash(head), khash(&rc_codes(tail)), khash(&rc_codes(head)), khash(tail)] { want.insert(key, Vec::new()); }
        }
        find_kmers(&region, region_start, w, &mut want);
        let mut round: Vec<(usize, bool, Vec<u8>)> = Vec::new();
        let mut left: Vec<Vec<u8>> = Vec::new();
        for t in rest {
            let r = rc_codes(&t);
            // Either end first: forward-strand anchors, then reverse-strand ones.
            let first = |s: &[u8]| want.get(&khash(s)).and_then(|v| v.first().copied());
            let hit = [&t, &r].into_iter().find_map(|s| first(&s[..w]).map(|p| (p, false, s.clone())))
                .or_else(|| [&t, &r].into_iter().find_map(|s| first(&rc_codes(&s[..w])).map(|p| (p, true, s.clone()))));
            match hit { Some(h) => round.push(h), None => left.push(t) }
        }
        rest = left;
        if round.is_empty() { break }
        round.sort_by_key(|x| x.0);
        region_start += region_len;
        region_len = round.iter().map(|(_, _, t)| t.len()).sum();
        rounds.push(round);
    }
    Plan { long, rounds, rest }
}

/// One tig as coded: its anchor (its first k - 1 bases are at that position
/// of the text, or their reverse complement if the flag is set), its length
/// (when encoding, the true one), and whether the decoder knows the length
/// (otherwise it is coded as end flags).
struct Slot { anchor: Option<(usize, bool)>, len: usize, known: bool }

/// Context of the end flag of an anchored tig before the base at `pos`: whether
/// a copy model has matched k - 1 bases exactly, the longest match, the position.
fn end_ctx(model: &Model, pos: usize) -> usize {
    let k = model.k;
    let e1 = model.models.iter().any(|cm| cm.m.active && cm.m.exact >= k - 1) as usize;
    let lb = model.longest_bucket();
    let pb = match pos - k { 0 => 0, 1..=2 => 1, 3..=15 => 2, 16..=63 => 3, 64..=127 => 4, 128..=253 => 5, _ => 6 };
    (e1 * 7 + lb) * 7 + pb
}

/// Codes (or, with `encode` false, decodes into `text`) the bases of the tigs
/// in `slots`. Fails if a base to encode is excluded (the tigs repeat a _k_-mer).
fn code_bases<C: Coder>(c: &mut C, model: &mut Model, text: &mut Vec<u8>, slots: &mut [Slot], encode: bool) -> Result<(), E> {
    let k = model.k;
    let w = k - 1;
    let mut end_tab = vec![C0; 2 * 7 * 7];
    let mut starts: Vec<usize> = Vec::with_capacity(slots.len() + 1);
    let mut i = 0usize;
    for slot in slots.iter_mut() {
        starts.push(i);
        let mut pos = 0usize;
        loop {
            if slot.known {
                if pos == slot.len { break }
            } else if pos >= k {
                if pos == 2 * k - 1 { slot.len = pos; break } // the longest short tig
                let x = end_ctx(model, pos);
                let bit = c.code((pos == slot.len) as u32, (end_tab[x] >> 16).clamp(1, 65535));
                update(&mut end_tab[x], bit);
                if bit == 1 { slot.len = pos; break }
            }
            let tig_of = |p: usize| starts.partition_point(|&x| x <= p) - 1;
            let (pr, ex) = model.predict(text, i, pos, &tig_of);
            let sym = match slot.anchor {
                Some((a, false)) if pos < w => text[a + pos],
                Some((a, true)) if pos < w => 3 - text[a + w - 1 - pos],
                _ => {
                    if encode && ex >> text[i] & 1 == 1 { return Err(format!("the tigs repeat a {}-mer", k).into()) }
                    code_base(c, pr, ex, if encode { text[i] } else { 0 })
                },
            };
            if !encode { text.push(sym); }
            model.update(text, i, pos, sym);
            i += 1;
            pos += 1;
        }
    }
    Ok(())
}

// Public interface -----------------------------------------------------------------

/// Writes `tigs` (A, C, G and T only, each at least `k` bases, longest first)
/// to `out`, see the module docs.
pub fn encode<W: Write>(
    tigs: &[Vec<u8>],
    k: usize,
    mut out: W,
) -> Result<W, E> {
    if k < 2 || k > u16::MAX as usize { return Err(format!("unsupported k = {}", k).into()) }
    if !tigs.windows(2).all(|w| w[0].len() >= w[1].len()) { return Err("tigs are not sorted longest first".into()) }
    if tigs.iter().any(|t| t.len() < k) { return Err(format!("a tig is shorter than k = {}", k).into()) }
    let mut codes: Vec<Vec<u8>> = Vec::with_capacity(tigs.len());
    for t in tigs {
        codes.push(t.iter().map(|&c| match c { b'A' => Ok(0), b'C' => Ok(1), b'G' => Ok(2), b'T' => Ok(3), _ => Err(format!("tigs contain {}", c as char)) }).collect::<Result<Vec<u8>, String>>()?);
    }
    let n: usize = codes.iter().map(|t| t.len()).sum();
    let plan = plan(codes, k);
    let mut model = Model::new(k, n);
    let mut enc = Encoder::new();
    let mut ll: Vec<usize> = plan.long.iter().map(|t| t.len()).collect();
    code_lengths(&mut enc, &mut ll, plan.long.len());
    let mut rl: Vec<usize> = plan.rest.iter().map(|t| t.len()).collect();
    code_lengths(&mut enc, &mut rl, plan.rest.len());
    let mut tpos = vec![C0; 128];
    let mut tstrand = C0;
    for round in &plan.rounds {
        let mut prev = 0usize;
        for (p, rev, _) in round {
            code_gamma(&mut enc, &mut tpos, (p - prev) as u64);
            prev = *p;
            let bit = enc.code(*rev as u32, (tstrand >> 16).clamp(1, 65535));
            update(&mut tstrand, bit);
        }
    }
    let mut text: Vec<u8> = Vec::with_capacity(n);
    let mut slots: Vec<Slot> = Vec::with_capacity(tigs.len());
    for t in &plan.long { text.extend_from_slice(t); slots.push(Slot { anchor: None, len: t.len(), known: true }); }
    for round in &plan.rounds { for (p, rev, t) in round { text.extend_from_slice(t); slots.push(Slot { anchor: Some((*p, *rev)), len: t.len(), known: false }); } }
    for t in &plan.rest { text.extend_from_slice(t); slots.push(Slot { anchor: None, len: t.len(), known: true }); }
    code_bases(&mut enc, &mut model, &mut text, &mut slots, true)?;
    let body = enc.finish();
    out.write_all(MAGIC)?;
    out.write_all(&(k as u16).to_le_bytes())?;
    out.write_all(&[OPT_END_FLAGS | OPT_RC_ANCHORS, plan.rounds.len() as u8])?;
    out.write_all(&(n as u64).to_le_bytes())?;
    for c in [plan.long.len(), plan.rest.len()].into_iter().chain(plan.rounds.iter().map(|r| r.len())) { out.write_all(&(c as u32).to_le_bytes())?; }
    out.write_all(&body)?;
    Ok(out)
}

/// Reads tigs written by [encode], in their coding order and orientation (the
/// same _k_-mers as the encoded tigs).
pub fn decode<R: Read>(
    mut input: R,
) -> Result<Vec<Vec<u8>>, E> {
    let mut data = Vec::new();
    input.read_to_end(&mut data)?;
    if data.len() < 24 || &data[..4] != MAGIC { return Err("not a tigs file (TAC5)".into()) }
    let k = u16::from_le_bytes([data[4], data[5]]) as usize;
    let (options, n_rounds) = (data[6], data[7] as usize);
    if options & OPT_END_FLAGS == 0 || options & !(OPT_END_FLAGS | OPT_RC_ANCHORS) != 0 { return Err("unsupported tigs file options".into()) }
    let rc_anchors = options & OPT_RC_ANCHORS != 0;
    let n = u64::from_le_bytes(data[8..16].try_into()?) as usize;
    let u32_at = |o: usize| -> Result<usize, E> { Ok(u32::from_le_bytes(data.get(o..o + 4).ok_or("truncated tigs file")?.try_into()?) as usize) };
    let (n_long, n_rest) = (u32_at(16)?, u32_at(20)?);
    let round_n: Vec<usize> = (0..n_rounds).map(|r| u32_at(24 + 4 * r)).collect::<Result<_, _>>()?;
    let mut model = Model::new(k, n);
    let mut dec = Decoder::new(&data[24 + 4 * n_rounds..]);
    let mut ll = Vec::new();
    code_lengths(&mut dec, &mut ll, n_long);
    let mut rl = Vec::new();
    code_lengths(&mut dec, &mut rl, n_rest);
    let mut tpos = vec![C0; 128];
    let mut tstrand = C0;
    let mut slots: Vec<Slot> = ll.iter().map(|&l| Slot { anchor: None, len: l, known: true }).collect();
    for &c in &round_n {
        let mut prev = 0usize;
        for _ in 0..c {
            let p = prev + code_gamma(&mut dec, &mut tpos, 0) as usize;
            prev = p;
            let rev = rc_anchors && {
                let bit = dec.code(0, (tstrand >> 16).clamp(1, 65535));
                update(&mut tstrand, bit);
                bit == 1
            };
            slots.push(Slot { anchor: Some((p, rev)), len: 0, known: false });
        }
    }
    slots.extend(rl.iter().map(|&l| Slot { anchor: None, len: l, known: true }));
    let mut text: Vec<u8> = Vec::with_capacity(n);
    code_bases(&mut dec, &mut model, &mut text, &mut slots, false)?;
    if text.len() != n { return Err("corrupt tigs file".into()) }
    let mut tigs = Vec::with_capacity(slots.len());
    let mut i = 0;
    for s in &slots {
        tigs.push(text[i..i + s.len].iter().map(|&c| b"ACGT"[c as usize]).collect());
        i += s.len;
    }
    Ok(tigs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::orientation::reverse_complement;

    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
        (0..n).map(|_| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b"ACGT"[(x >> 33) as usize % 4] }).collect()
    }

    fn canonical_sorted(tigs: &[Vec<u8>]) -> Vec<Vec<u8>> {
        let mut out: Vec<Vec<u8>> = tigs.iter().map(|t| { let r = reverse_complement(t); if r < *t { r } else { t.clone() } }).collect();
        out.sort();
        out
    }

    /// `s` followed by `tail`, whose first base is changed if it equals `avoid`.
    fn branch(s: &[u8], mut tail: Vec<u8>, avoid: u8) -> Vec<u8> {
        if tail[0] == avoid { tail[0] = if avoid == b'A' { b'C' } else { b'A' }; }
        let mut out = s.to_vec();
        out.extend_from_slice(&tail);
        out
    }

    #[test]
    fn roundtrip_with_anchors() {
        let k = 31;
        // A spectrum-preserving set: two long tigs; a SNP variant of one, a
        // branch off the other (written on its reverse strand), and a tig whose
        // only anchor is on the reverse strand of the other (anchored in the
        // long tigs); a branch off the SNP variant (anchored in it); unrelated
        // tigs.
        let a = random_seq(900, 1);
        let b = random_seq(700, 2);
        let mut snp = a[300..300 + 2 * k - 1].to_vec();
        snp[k - 1] = if snp[k - 1] == b'A' { b'C' } else { b'A' };
        let branch_rc = reverse_complement(&branch(&b[100..100 + k - 1], random_seq(15, 3), b[100 + k - 1]));
        let reverse_only = branch(&reverse_complement(&b[400..400 + k - 1]), random_seq(10, 7), reverse_complement(&b[399..400])[0]);
        let second = branch(&snp[1..k], random_seq(12, 4), snp[k]);
        let mut tigs = vec![a, b, snp, branch_rc, reverse_only, second, random_seq(2 * k - 5, 5), random_seq(k, 6)];
        tigs.sort_by_key(|t| std::cmp::Reverse(t.len()));
        let codes: Vec<Vec<u8>> = tigs.iter().map(|t| t.iter().map(|&c| b"ACGT".iter().position(|&x| x == c).unwrap() as u8).collect()).collect();
        let p = plan(codes, k);
        // Four anchored tigs in two rounds, one of them on the reverse strand.
        assert_eq!(p.rounds.iter().map(|r| r.len()).collect::<Vec<_>>(), vec![3, 1]);
        assert_eq!(p.rounds[0].iter().filter(|x| x.1).count(), 1);
        let bytes = encode(&tigs, k, Vec::new()).unwrap();
        assert_eq!((bytes[6], bytes[7]), (OPT_END_FLAGS | OPT_RC_ANCHORS, 2));
        let decoded = decode(bytes.as_slice()).unwrap();
        assert_eq!(canonical_sorted(&decoded), canonical_sorted(&tigs));
    }

    #[test]
    fn rejects_repeated_kmers() {
        let k = 31;
        let a = random_seq(200, 7);
        assert!(encode(&[a.clone(), a[10..80].to_vec()], k, Vec::new()).is_err());
    }

    #[test]
    fn rejects_bad_input() {
        assert!(encode(&[b"ACGT".to_vec(), b"ACGTACGT".to_vec()], 3, Vec::new()).is_err());
        assert!(encode(&[b"ACGN".to_vec()], 3, Vec::new()).is_err());
        assert!(encode(&[b"AC".to_vec()], 3, Vec::new()).is_err());
        assert!(decode(&b"TAC4"[..]).is_err());
    }
}
