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
//! Index tigs: the _k_-mer set of the compression index, chosen from the
//! input assemblies and written as compressed simplitigs.
//!
//! The index is shipped as its tigs because building the SBWT from them is
//! deterministic. [build] chooses the _k_-mers in four steps:
//!
//! 1. **Core**: the canonical _k_-mers present in at least
//!    [TigOpts::min_assemblies] of the inputs. They are counted exactly with
//!    canonical 128-bit rolling hashes, one hash partition at a time to bound
//!    the memory, and an SBWT is built from them.
//! 2. **Gap filling**: every input is streamed through the core SBWT. A gap
//!    is a maximal run of _g_ windows not in the core with a core _k_-mer on
//!    both sides; a run at a contig end is a tip. A gap is a *join* if the
//!    _k_-mer before it has no successor and the one after it no predecessor
//!    in the core (both are tig ends), *one-sided* if exactly one of them is
//!    a dead end, and a *bubble* otherwise. The _k_-mers of the gaps and
//!    tips whose _g_ is within the ranges in [GapRules] are added. With
//!    [GapRules::drop_repairable_bubbles], a bubble the archive can store as
//!    edits is left out: one whose window (the gap with the core _k_-mers on
//!    both sides) [call_snps](crate::graph::snps::call_snps) edits into core
//!    _k_-mers only.
//! 3. **Simplitigs**: an SBWT is built from the core and the gap _k_-mers
//!    (this is the final index) and its de Bruijn graph is covered by greedy
//!    simplitigs. The inputs are scanned in order; each tig starts at the
//!    first _k_-mer not used yet and is extended to the right, then to the
//!    left, through unused _k_-mers. While a tig still follows the input it
//!    was seeded from, it takes that input's next base if it can; otherwise
//!    it tries A, C, G and T in that order.
//! 4. **Layout**: the tigs are put on the strand [StrandChooser] picks for
//!    them in the order they were made, and sorted longest first. This
//!    changes only their compressed size, not their _k_-mers: the long,
//!    core-like tigs come first and the short variant tigs are coded against
//!    them.
//!
//! [write_tigs] stores the tigs with a DNA context-mixing coder (see
//! [crate::coder::tigs]), [read_tigs] reads them back (in the coder's order and
//! orientation: the same _k_-mers).
//!
use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{Read, Write};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU32, Ordering};

use rayon::prelude::*;
use crate::walk::{IndexWalk, LcsIndex};
use sbwt::sbwt_index_variant::SbwtIndexVariant;
use sbwt::{LcsArray, SbwtIndex, StreamingIndex, SubsetMatrix};

use crate::graph::orientation::{reverse_complement, StrandChooser};

use crate::coder;

type E = Box<dyn std::error::Error>;

/// Inclusive range of gap lengths (in _k_-mer windows) whose _k_-mers are
/// added; `None` adds none.
pub type GapRange = Option<(usize, usize)>;

/// Which gaps against the core are filled, by class.
#[derive(Clone, Debug)]
pub struct GapRules {
    /// Both flanks are tig ends.
    pub join: GapRange,
    /// Exactly one flank is a tig end.
    pub one_sided: GapRange,
    /// Neither flank is a tig end.
    pub bubble: GapRange,
    /// Runs at a contig end (their valid windows are added).
    pub tip: GapRange,
    /// Leave out the bubbles that SNP and indel editing repairs against the
    /// core (see the module docs).
    pub drop_repairable_bubbles: bool,
}

impl Default for GapRules {
    /// All joins and one-sided gaps; bubbles of at least 350 windows that
    /// editing does not repair; tips of at least 250. The best rules found at
    /// 1,991 _E. coli_ assemblies with a core of _k_-mers in at least 5 of
    /// them, for archives with SNP and indel editing (INDEX_JOURNAL.md).
    fn default() -> Self {
        GapRules {
            join: Some((1, usize::MAX)),
            one_sided: Some((1, usize::MAX)),
            bubble: Some((350, usize::MAX)),
            tip: Some((250, usize::MAX)),
            drop_repairable_bubbles: true,
        }
    }
}

fn in_range(range: &GapRange, g: usize) -> bool {
    range.is_some_and(|(lo, hi)| g >= lo && g <= hi)
}

/// Parameters for [build].
#[derive(Clone, Debug)]
pub struct TigOpts {
    /// Core: canonical _k_-mers in at least this many inputs.
    pub min_assemblies: u32,
    /// Number of hash partitions the core is counted in (more partitions use
    /// less memory but read the inputs more times).
    pub partitions: u64,
    /// Gaps against the core whose _k_-mers are added.
    pub gaps: GapRules,
    /// Sampling width of [StrandChooser] for orienting the tigs.
    pub orient_sampling: u64,
    /// SBWT construction; `k`, `num_threads`, `mem_gb`, `temp_dir` and
    /// `prefix_precalc` are used, reverse complements and select support are
    /// always added.
    pub build_opts: kbo::BuildOpts,
}

impl Default for TigOpts {
    fn default() -> Self {
        let mut build_opts = kbo::BuildOpts::default();
        build_opts.k = 255;
        TigOpts {
            min_assemblies: 5,
            partitions: 16,
            gaps: GapRules::default(),
            orient_sampling: 8,
            build_opts,
        }
    }
}

/// Counts from [build].
#[derive(Clone, Debug, Default)]
pub struct TigStats {
    /// Canonical _k_-mers in the core.
    pub core_kmers: usize,
    /// Gap occurrences whose _k_-mers were added, by class: join, one-sided,
    /// bubble, tip.
    pub gaps: [usize; 4],
    /// Windows in those gap occurrences, by class.
    pub gap_windows: [usize; 4],
    /// Bubble occurrences in range left out as repairable.
    pub repairable_bubbles: usize,
    /// Canonical _k_-mers in the final index.
    pub kmers: usize,
    /// Tigs reversed by the layout.
    pub reversed: usize,
}

/// Result of [build]: the tigs and the index built from them.
pub struct Tigs {
    /// The tigs, oriented and longest first.
    pub tigs: Vec<Vec<u8>>,
    /// SBWT of the tigs (with reverse complements); the same as building it
    /// from [Tigs::tigs].
    pub sbwt: SbwtIndexVariant,
    /// LCS array of [Tigs::sbwt].
    pub lcs: LcsArray,
    pub stats: TigStats,
}

// Canonical rolling hashes ---------------------------------------------------

/// The two 64-bit halves of a hash, combined. Both halves have the same
/// lowest bit (the code of the last base), so one is rotated first.
fn fold(x: u128) -> u64 {
    x as u64 ^ ((x >> 64) as u64).rotate_left(1)
}

/// The keys are already hashes: use them folded to 64 bits.
#[derive(Default)]
struct Fold(u64);
impl Hasher for Fold {
    fn finish(&self) -> u64 { self.0 }
    fn write(&mut self, bytes: &[u8]) { for &b in bytes { self.0 = (self.0 << 8) | b as u64; } }
    fn write_u128(&mut self, x: u128) { self.0 = fold(x); }
}
type H = BuildHasherDefault<Fold>;

const B1: u64 = 0x9E3779B97F4A7C15;
const B2: u64 = 0xC2B2AE3D27D4EB4F;

fn inv(b: u64) -> u64 { let mut x = b; for _ in 0..6 { x = x.wrapping_mul(2u64.wrapping_sub(b.wrapping_mul(x))); } x }
fn pow(b: u64, e: usize) -> u64 { (0..e).fold(1u64, |a, _| a.wrapping_mul(b)) }
fn code(c: u8) -> Option<u64> { match c { b'A' => Some(0), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None } }
fn is_acgt(c: u8) -> bool { code(c).is_some() }

/// Canonical hash of every _k_-mer window of `seq` (None for windows with a
/// character other than A, C, G or T). Two 64-bit polynomial hashes per
/// strand; the canonical hash is the smaller strand's.
fn hashes(seq: &[u8], k: usize) -> Vec<Option<u128>> {
    let (p1, p2) = (pow(B1, k - 1), pow(B2, k - 1));
    let (i1, i2) = (inv(B1), inv(B2));
    let mut out = Vec::with_capacity(seq.len().saturating_sub(k - 1));
    let (mut f1, mut f2, mut r1, mut r2) = (0u64, 0u64, 0u64, 0u64);
    let mut valid = 0usize; // consecutive ACGT bases ending here
    for i in 0..seq.len() {
        match code(seq[i]) {
            None => { valid = 0; f1 = 0; f2 = 0; r1 = 0; r2 = 0; },
            Some(c) => {
                if valid >= k {
                    let o = code(seq[i - k]).unwrap();
                    f1 = f1.wrapping_sub(o.wrapping_mul(p1)); f2 = f2.wrapping_sub(o.wrapping_mul(p2));
                    r1 = r1.wrapping_sub(3 - o).wrapping_mul(i1); r2 = r2.wrapping_sub(3 - o).wrapping_mul(i2);
                    f1 = f1.wrapping_mul(B1).wrapping_add(c); f2 = f2.wrapping_mul(B2).wrapping_add(c);
                    r1 = r1.wrapping_add((3 - c).wrapping_mul(p1)); r2 = r2.wrapping_add((3 - c).wrapping_mul(p2));
                } else {
                    f1 = f1.wrapping_mul(B1).wrapping_add(c); f2 = f2.wrapping_mul(B2).wrapping_add(c);
                    r1 = r1.wrapping_add((3 - c).wrapping_mul(pow(B1, valid))); r2 = r2.wrapping_add((3 - c).wrapping_mul(pow(B2, valid)));
                }
                valid += 1;
            },
        }
        if i + 1 >= k {
            out.push(if valid >= k {
                let f = ((f1 as u128) << 64) | f2 as u128;
                let r = ((r1 as u128) << 64) | r2 as u128;
                Some(f.min(r))
            } else { None });
        }
    }
    out
}

// Steps ----------------------------------------------------------------------

/// Canonical _k_-mers present in at least `t` of the inputs.
fn count_core<T: Sync>(
    inputs: &[T],
    read: &(dyn Fn(&T) -> Vec<Vec<u8>> + Sync),
    k: usize,
    t: u32,
    partitions: u64,
) -> HashSet<u128, H> {
    let mut core: HashSet<u128, H> = HashSet::default();
    for p in 0..partitions {
        let shards: Vec<Mutex<HashMap<u128, u32, H>>> = (0..256).map(|_| Mutex::new(HashMap::default())).collect();
        inputs.par_iter().for_each(|input| {
            let mut hs: Vec<u128> = Vec::new();
            for seq in read(input) {
                hs.extend(hashes(&seq, k).into_iter().flatten().filter(|&h| fold(h) % partitions == p));
            }
            hs.sort_unstable();
            hs.dedup();
            for group in hs.chunk_by(|a, b| (a >> 120) == (b >> 120)) {
                let mut counts = shards[(group[0] >> 120) as usize].lock().unwrap();
                for &h in group { *counts.entry(h).or_insert(0) += 1; }
            }
        });
        for shard in shards {
            core.extend(shard.into_inner().unwrap().into_iter().filter(|&(_, m)| m >= t).map(|(h, _)| h));
        }
    }
    core
}

/// The _k_-mers of `set` as runs of the inputs, each _k_-mer once. Empties
/// `set`.
fn kmer_runs<T: Sync>(
    inputs: &[T],
    read: &(dyn Fn(&T) -> Vec<Vec<u8>> + Sync),
    k: usize,
    set: &mut HashSet<u128, H>,
) -> Vec<Vec<u8>> {
    let mut runs: Vec<Vec<u8>> = Vec::new();
    let chunk = rayon::current_num_threads().max(1) * 2;
    for batch in inputs.chunks(chunk) {
        if set.is_empty() { break }
        let hashed: Vec<Vec<(Vec<u8>, Vec<Option<u128>>)>> = batch.par_iter()
            .map(|input| read(input).into_iter().map(|seq| { let hs = hashes(&seq, k); (seq, hs) }).collect())
            .collect();
        for (seq, hs) in hashed.into_iter().flatten() {
            let mut run: Option<usize> = None;
            for (w, h) in hs.iter().enumerate().chain(std::iter::once((hs.len(), &None))) {
                if h.is_some_and(|h| set.remove(&h)) {
                    if run.is_none() { run = Some(w); }
                } else if let Some(s) = run.take() {
                    runs.push(seq[s..w + k - 1].to_vec());
                }
            }
        }
    }
    runs
}

/// SBWT with reverse complements and select support.
fn build_index(
    seqs: &[Vec<u8>],
    opts: &kbo::BuildOpts,
) -> (SbwtIndexVariant, LcsArray) {
    let mut opts = opts.clone();
    opts.add_revcomp = true;
    opts.build_select = true;
    kbo::build(seqs, opts)
}

/// Colex rank of every _k_-mer window of `seq` (None if not in the index).
fn colex_ranks(
    index: &StreamingIndex<SbwtIndex<SubsetMatrix>, LcsArray>,
    seq: &[u8],
    k: usize,
) -> Vec<Option<usize>> {
    if seq.len() < k { return Vec::new() }
    index.matching_statistics(seq)[k - 1..].iter()
        .map(|(len, range)| if *len == k { Some(range.start) } else { None })
        .collect()
}

/// Whether editing (SNPs and small indels, see
/// [call_snps](crate::graph::snps::call_snps)) turns `window`, which starts
/// and ends with a core _k_-mer, into core _k_-mers only.
fn repairable(
    window: &[u8],
    core: &SbwtIndexVariant,
    lcs: &LcsArray,
    index: &StreamingIndex<SbwtIndex<SubsetMatrix>, LcsArray>,
    k: usize,
) -> bool {
    let edits: Vec<crate::graph::snps::Edit> = crate::graph::snps::call_snps(window, core, lcs).into_iter().map(|called| called.edit).collect();
    if edits.is_empty() {
        return false
    }
    crate::graph::snps::apply_edits(window, edits)
        .is_some_and(|(edited, _)| edited.len() >= k && index.matching_statistics(&edited)[k - 1..].iter().all(|m| m.0 == k))
}

/// The windows of the gaps of one input that `rules` fill, as sequences.
/// `core` is the core index (for [GapRules::drop_repairable_bubbles]);
/// `counts` holds the occurrences and windows per class, then the bubbles
/// left out as repairable.
#[allow(clippy::too_many_arguments)]
fn input_gaps(
    seqs: &[Vec<u8>],
    index: &StreamingIndex<SbwtIndex<SubsetMatrix>, LcsArray>,
    walk: &LcsIndex,
    core: (&SbwtIndexVariant, &LcsArray),
    k: usize,
    rules: &GapRules,
    gaps: &mut Vec<Vec<u8>>,
    counts: &mut [usize; 9],
) {
    for seq in seqs {
        let ranks = colex_ranks(index, seq, k);
        let nw = ranks.len();
        let mut w = 0;
        while w < nw {
            if ranks[w].is_some() { w += 1; continue }
            let w0 = w;
            while w < nw && ranks[w].is_none() { w += 1; }
            let w1 = w - 1;
            let g = w1 - w0 + 1;
            let span = &seq[w0..w1 + k];
            if w0 == 0 || w1 == nw - 1 {
                if in_range(&rules.tip, g) {
                    // the valid windows of the run
                    for piece in span.split(|&c| !is_acgt(c)).filter(|p| p.len() >= k) {
                        gaps.push(piece.to_vec());
                    }
                    counts[3] += 1; counts[7] += g;
                }
                continue
            }
            if span.iter().any(|&c| !is_acgt(c)) { continue }
            let (left, right) = (ranks[w0 - 1].unwrap(), ranks[w1 + 1].unwrap());
            let dead = (walk.out_degree(left) == 0) as usize + (walk.in_degree(right) == 0) as usize;
            let class = match dead { 2 => 0, 1 => 1, _ => 2 };
            let range = [&rules.join, &rules.one_sided, &rules.bubble][class];
            if in_range(range, g) {
                if class == 2 && rules.drop_repairable_bubbles && repairable(&seq[w0 - 1..w1 + k + 1], core.0, core.1, index, k) {
                    counts[8] += 1;
                    continue
                }
                gaps.push(span.to_vec());
                counts[class] += 1; counts[4 + class] += g;
            }
        }
    }
}

/// Greedy simplitigs covering the index, seeded in input order.
fn simplitigs<T: Sync>(
    inputs: &[T],
    read: &(dyn Fn(&T) -> Vec<Vec<u8>> + Sync),
    sbwt: &SbwtIndex<SubsetMatrix>,
    lcs: &LcsArray,
    runs: &[Vec<u8>],
) -> Vec<Vec<u8>> {
    let k = sbwt.k();
    let n = sbwt.n_sets();
    let index = StreamingIndex::new(sbwt, lcs);
    let walk = LcsIndex { sbwt, lcs };

    // Colex rank of the reverse complement of every k-mer.
    let rc: Vec<AtomicU32> = (0..n).map(|_| AtomicU32::new(u32::MAX)).collect();
    runs.par_iter().for_each(|run| {
        let f = colex_ranks(&index, run, k);
        let b = colex_ranks(&index, &reverse_complement(run), k);
        let m = f.len();
        for j in 0..m {
            let (cf, cb) = (f[j].unwrap(), b[m - 1 - j].unwrap());
            rc[cf].store(cb as u32, Ordering::Relaxed);
            rc[cb].store(cf as u32, Ordering::Relaxed);
        }
    });
    let rc: Vec<u32> = rc.into_iter().map(|x| x.into_inner()).collect();

    let mut used = vec![0u64; n.div_ceil(64)];
    let is_used = |used: &[u64], v: usize| used[v / 64] >> (v % 64) & 1 == 1;
    let set_used = |used: &mut [u64], v: usize| used[v / 64] |= 1 << (v % 64);
    let complement = |c: u8| match c { b'A' => b'T', b'C' => b'G', b'G' => b'C', _ => b'A' };
    let mut neighbors: Vec<(usize, u8)> = Vec::new();
    // Extends from `start` through unused k-mers; returns the labels taken.
    // While the tig follows its seed input, `seed(j)` (the input's label at
    // step j) is tried first; otherwise the labels in `order`.
    let mut extend = |used: &mut Vec<u64>, start: usize, order: &[u8; 4], seed: &dyn Fn(usize) -> Option<u8>| -> Vec<u8> {
        let mut labels = Vec::new();
        let mut cur = start;
        let mut on_seed = true;
        loop {
            neighbors.clear();
            walk.push_successors(cur, &mut neighbors);
            let s = if on_seed { seed(labels.len()) } else { None };
            let next = s.iter().chain(order.iter().filter(|&&c| Some(c) != s))
                .find_map(|&c| neighbors.iter().find(|(v, l)| *l == c && !is_used(used, *v)).map(|(v, l)| (*v, *l)));
            match next {
                Some((v, label)) => {
                    set_used(used, v);
                    set_used(used, rc[v] as usize);
                    labels.push(label);
                    on_seed &= Some(label) == s;
                    cur = v;
                },
                None => break,
            }
        }
        labels
    };

    let mut tigs: Vec<Vec<u8>> = Vec::new();
    let chunk = rayon::current_num_threads().max(1) * 2;
    for batch in inputs.chunks(chunk) {
        let ranked: Vec<Vec<(Vec<u8>, Vec<Option<usize>>)>> = batch.par_iter()
            .map(|input| read(input).into_iter().map(|seq| { let r = colex_ranks(&index, &seq, k); (seq, r) }).collect())
            .collect();
        for (seq, ranks) in ranked.into_iter().flatten() {
            for (w, v) in ranks.into_iter().enumerate() {
                let Some(v) = v else { continue };
                if is_used(&used, v) { continue }
                set_used(&mut used, v);
                set_used(&mut used, rc[v] as usize);
                // Right, then left as the right extension of the reverse
                // complement: prepending A, C, G, T is appending T, G, C, A.
                let right = extend(&mut used, v, b"ACGT", &|j| seq.get(w + k + j).copied());
                let left = extend(&mut used, rc[v] as usize, b"TGCA", &|j| (j < w).then(|| seq[w - 1 - j]).filter(|c| b"ACGT".contains(c)).map(complement));
                let mut tig: Vec<u8> = left.iter().rev().map(|&c| complement(c)).collect();
                tig.extend_from_slice(&seq[w..w + k]);
                tig.extend_from_slice(&right);
                tigs.push(tig);
            }
        }
    }
    tigs
}

/// Puts each tig on the strand `StrandChooser` picks in the given order and
/// sorts them longest first (stable). Returns the number reversed.
fn layout(
    tigs: &mut [Vec<u8>],
    sampling: u64,
) -> usize {
    let mut chooser = StrandChooser::new(sampling);
    let mut reversed = 0;
    for tig in tigs.iter_mut() {
        if chooser.choose(tig) {
            *tig = reverse_complement(tig);
            reversed += 1;
        }
    }
    tigs.sort_by_key(|t| std::cmp::Reverse(t.len()));
    reversed
}

/// Chooses the index _k_-mers from `inputs` and returns them as tigs with
/// the SBWT built from them, see the [module docs](self). `read` returns the
/// sequences of one input; `progress` is called after each step.
pub fn build<T: Sync>(
    inputs: &[T],
    read: &(dyn Fn(&T) -> Vec<Vec<u8>> + Sync),
    opts: &TigOpts,
    progress: &dyn Fn(&str),
) -> Tigs {
    let k = opts.build_opts.k;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(opts.build_opts.num_threads.max(1)).build().unwrap();
    let mut stats = TigStats::default();

    let mut core = pool.install(|| count_core(inputs, read, k, opts.min_assemblies, opts.partitions));
    stats.core_kmers = core.len();
    progress(&format!("core: {} canonical k-mers in at least {} inputs", core.len(), opts.min_assemblies));
    let mut runs = pool.install(|| kmer_runs(inputs, read, k, &mut core));
    drop(core);
    let (core_sbwt, core_lcs) = build_index(&runs, &opts.build_opts);
    progress("core: SBWT built");

    let mut counts = [0usize; 9];
    let gaps: Vec<Vec<u8>> = pool.install(|| {
        let SbwtIndexVariant::SubsetMatrix(core) = &core_sbwt;
        let index = StreamingIndex::new(core, &core_lcs);
        let walk = LcsIndex { sbwt: core, lcs: &core_lcs };
        let per_input: Vec<(Vec<Vec<u8>>, [usize; 9])> = inputs.par_iter().map(|input| {
            let (mut gaps, mut counts) = (Vec::new(), [0usize; 9]);
            input_gaps(&read(input), &index, &walk, (&core_sbwt, &core_lcs), k, &opts.gaps, &mut gaps, &mut counts);
            (gaps, counts)
        }).collect();
        per_input.into_iter().flat_map(|(gaps, c)| { for i in 0..9 { counts[i] += c[i]; } gaps }).collect()
    });
    drop(core_sbwt);
    drop(core_lcs);
    stats.gaps.copy_from_slice(&counts[0..4]);
    stats.gap_windows.copy_from_slice(&counts[4..8]);
    stats.repairable_bubbles = counts[8];
    progress(&format!("gaps filled (join, one-sided, bubble, tip): {:?}, windows {:?}; repairable bubbles left out: {}", stats.gaps, stats.gap_windows, stats.repairable_bubbles));

    runs.extend(gaps);
    let (sbwt, lcs) = build_index(&runs, &opts.build_opts);
    progress("index: SBWT built");
    let SbwtIndexVariant::SubsetMatrix(index_sbwt) = &sbwt;
    stats.kmers = index_sbwt.n_kmers() / 2;
    let mut tigs = pool.install(|| simplitigs(inputs, read, index_sbwt, &lcs, &runs));
    drop(runs);
    progress(&format!("{} simplitigs", tigs.len()));
    stats.reversed = layout(&mut tigs, opts.orient_sampling);

    Tigs { tigs, sbwt, lcs, stats }
}

// Tig files -------------------------------------------------------------------

/// Writes `tigs` (_k_-mers of length `k`, each tig at least `k` bases,
/// longest first, no _k_-mer twice) with the DNA coder ([coder::tigs::encode]).
pub fn write_tigs<W: Write>(
    tigs: &[Vec<u8>],
    k: usize,
    out: W,
) -> Result<W, E> {
    coder::tigs::encode(tigs, k, out)
}

/// Reads tigs written by [write_tigs], in the coder's order and orientation:
/// the same _k_-mers, so the same SBWT.
pub fn read_tigs<R: Read>(
    input: R,
) -> Result<Vec<Vec<u8>>, E> {
    coder::tigs::decode(input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hashes_are_canonical() {
        let s = b"ACGTTGCAAGGCTTACGATCGATCGGATC";
        let k = 7;
        let f = hashes(s, k);
        let r = hashes(&reverse_complement(s), k);
        let mut rr = r.clone();
        rr.reverse();
        assert_eq!(f, rr);
        assert!(hashes(b"ACGTNACGTACG", 5).iter().take(5).all(|h| h.is_none()));
    }

    #[test]
    fn gap_ranges() {
        let rules = GapRules::default();
        assert!(in_range(&rules.join, 1));
        assert!(in_range(&rules.one_sided, 1));
        assert!(!in_range(&rules.bubble, 349));
        assert!(in_range(&rules.bubble, 350));
        assert!(!in_range(&rules.tip, 249));
        assert!(in_range(&rules.tip, 250));
        assert!(rules.drop_repairable_bubbles);
        assert!(!in_range(&None, 5));
    }

    // Random sequence for the end-to-end test.
    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
        (0..n).map(|_| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b"ACGT"[(x >> 33) as usize % 4] }).collect()
    }

    /// Canonical k-mers of a set of sequences.
    fn kmer_set(seqs: &[Vec<u8>], k: usize) -> HashSet<Vec<u8>> {
        let mut out = HashSet::new();
        for s in seqs {
            for w in s.windows(k) {
                let r = reverse_complement(w);
                out.insert(if r.as_slice() < w { r } else { w.to_vec() });
            }
        }
        out
    }

    #[test]
    fn repairable_bubbles_are_left_out() {
        let k = 31;
        // Two inputs are the backbone (the core); a third has a SNP (a bubble
        // the edits repair) and a fourth an insertion of 40 bases (a bubble
        // they do not).
        let backbone = random_seq(2000, 4);
        let mut with_snp = backbone.clone();
        with_snp[700] = if backbone[700] == b'A' { b'C' } else { b'A' };
        let mut with_insert = backbone[..1300].to_vec();
        with_insert.extend_from_slice(&random_seq(40, 5));
        with_insert.extend_from_slice(&backbone[1300..]);
        let inputs: Vec<Vec<Vec<u8>>> = vec![vec![backbone.clone()], vec![backbone.clone()], vec![with_snp.clone()], vec![with_insert.clone()]];
        let read = |x: &Vec<Vec<u8>>| x.clone();
        let mut opts = TigOpts::default();
        opts.min_assemblies = 2;
        opts.partitions = 2;
        opts.build_opts.k = k;
        opts.build_opts.num_threads = 2;
        opts.gaps = GapRules { join: Some((1, usize::MAX)), one_sided: Some((1, usize::MAX)), bubble: Some((1, usize::MAX)), tip: Some((1, usize::MAX)), drop_repairable_bubbles: false };
        let all = build(&inputs, &read, &opts, &|_| {});
        assert_eq!((all.stats.gaps[2], all.stats.repairable_bubbles), (2, 0));
        opts.gaps.drop_repairable_bubbles = true;
        let dropped = build(&inputs, &read, &opts, &|_| {});
        assert_eq!((dropped.stats.gaps[2], dropped.stats.repairable_bubbles), (1, 1));
        // The SNP's k-mers are left out, the insertion's kept.
        let (snp_kmers, insert_kmers) = (kmer_set(&[with_snp], k), kmer_set(&[with_insert], k));
        let index = kmer_set(&dropped.tigs, k);
        assert!(snp_kmers.difference(&kmer_set(&[backbone.clone()], k)).all(|x| !index.contains(x)));
        assert!(insert_kmers.is_subset(&index));
        assert!(snp_kmers.is_subset(&kmer_set(&all.tigs, k)));
    }

    #[test]
    fn build_core_gaps_and_layout() {
        let k = 15;
        // Three inputs share a backbone; one carries an insertion (a bubble of
        // its own k-mers) and another a different tail (a tip).
        let backbone = random_seq(400, 1);
        let insert = random_seq(60, 2);
        let tail = random_seq(80, 3);
        let mut with_insert = backbone[..200].to_vec();
        with_insert.extend_from_slice(&insert);
        with_insert.extend_from_slice(&backbone[200..]);
        let mut with_tail = backbone.clone();
        with_tail.extend_from_slice(&tail);
        let inputs: Vec<Vec<Vec<u8>>> = vec![vec![backbone.clone()], vec![with_insert.clone()], vec![reverse_complement(&with_tail)]];
        let read = |x: &Vec<Vec<u8>>| x.clone();

        let mut opts = TigOpts::default();
        opts.min_assemblies = 3;
        opts.partitions = 3;
        opts.build_opts.k = k;
        opts.build_opts.num_threads = 2;
        opts.gaps = GapRules { join: Some((1, usize::MAX)), one_sided: Some((1, usize::MAX)), bubble: Some((1, usize::MAX)), tip: Some((1, usize::MAX)), drop_repairable_bubbles: false };
        let result = build(&inputs, &read, &opts, &|_| {});

        // Core = the backbone k-mers outside the insertion breakpoint (all three
        // inputs have them); every input k-mer ends up in the index.
        let all = kmer_set(&[backbone.clone(), with_insert.clone(), with_tail.clone()], k);
        assert_eq!(kmer_set(&result.tigs, k), all);
        assert_eq!(result.stats.kmers, all.len());
        // Each k-mer once: total windows equals the k-mer count.
        let windows: usize = result.tigs.iter().map(|t| t.len() - k + 1).sum();
        assert_eq!(windows, all.len());
        // Longest first.
        assert!(result.tigs.windows(2).all(|w| w[0].len() >= w[1].len()));
        // The k-mers across the insertion point are in two inputs only, so
        // every input has a gap there whose flanks are both core tig ends
        // (three joins); the tail is a tip.
        assert_eq!(result.stats.gaps, [3, 0, 0, 1]);
        // The tig file restores the same k-mers.
        let bytes = write_tigs(&result.tigs, k, Vec::new()).unwrap();
        assert_eq!(kmer_set(&read_tigs(bytes.as_slice()).unwrap(), k), all);
    }
}
