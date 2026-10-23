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
//! Context-mixing model for the choice section.
//!
//! The choices of all colours are coded as one stream with the arithmetic
//! coder ([crate::codecs::arith]). Colours are coded in order, each by walking
//! its path through the graph, so the decoder needs the graph first. At every
//! node with `n > 1` outgoing edges the rank of the taken edge is coded as
//! binary decisions "is it candidate _c_?" over the candidates in order:
//! 1. the predictions of the colour's reference colours (up to [N_REFS]
//!    earlier colours, by score: their choice at this node);
//! 2. the choice of the most recent colour at this node;
//! 3. the choice of the previous colour (colour - 1) at this node;
//! 4. the ranks not yet listed in order, for `n` ≤ [MAX_LISTED]; for larger
//!    `n` the rank through an adaptive binary tree of its bits instead.
//!
//! A candidate that is the only one left is implied. The references start
//! from the most recent colour at the first node with a history. A reference
//! that predicts right scores +1, one that predicts wrong has its score
//! halved; if none was right, the most recent colour that made the same
//! choice joins (replacing the lowest score).
//!
//! The probability of every decision mixes, in the logistic domain:
//! - an adaptive counter per edge (taken when offered);
//! - for copy candidates, a counter per candidate kind and position;
//! - a counter by the number of references predicting the candidate;
//! - the share of the earlier visits to the node that took the edge;
//! - a copy flag and a bias.
//!
//! Three mixers select their weight sets by the candidate class (top
//! reference, other reference, other copy, plain rank) and: the out-degree;
//! the votes and the top reference's score; the visits to the node and the
//! candidate's position. A final mixer per class combines them, and an APM
//! ([Apm], fixed rate) per class refines the result. Everything is integer arithmetic, so the
//! decoder repeats the predictions exactly.
//!
//! The format is the arithmetic-coded stream alone: the number of colours
//! comes from the file header, and the decoder stops at every colour's sink.
use std::collections::HashMap;

use crate::codecs::arith::{code_tree, update, Coder, Decoder, Encoder, C0};
use crate::codecs::mix::{Apm, ApmRate, Tables, SMAX};

use super::choices::{ChoiceStream, ChoiceWriter, Choices};
use super::{GraphWalk, StoredGraph, SINK_NODE, SOURCE_NODE};

/// Reference colours per colour.
const N_REFS: usize = 8;
/// Largest out-degree whose ranks are listed as candidates.
const MAX_LISTED: u32 = 16;
/// Recent colours kept per edge.
const N_TOOK: usize = 4;
/// Candidate kinds: the references, the most recent colour, the previous
/// colour, a plain rank.
const KIND_LAST: usize = N_REFS;
const KIND_PREV: usize = N_REFS + 1;
const KIND_RANK: usize = N_REFS + 2;
const N_KINDS: usize = N_REFS + 3;
/// Candidate classes: top reference, other reference, other copy, rank.
const N_CLASSES: usize = 4;

/// Mixer inputs and weight sets.
const N_IN: usize = 6;
const SETS_A: usize = N_CLASSES * 3;
const SETS_B: usize = N_CLASSES * 8 * 8;
const SETS_C: usize = N_CLASSES * 12 * 4;
const N_FIN: usize = 4;
/// Mixer weight 1.0; learning rates as dw = x * err * LR >> 16 (0.02 and
/// 0.01 with x in 1/256 nat and err in 1/65536).
const W_ONE: i64 = 1 << 20;
const LR: i64 = 82;
const LR_FINAL: i64 = 41;
/// Bias input (0.3) and copy flag (1.0) in 1/256 nat.
const BIAS: i32 = 77;
const FLAG: i32 = 256;
/// A reference's score of 1.
const SCORE_ONE: u32 = 256;
/// The APM's fixed update rate (0.02 in 1/65536).
const APM_RATE: i64 = 1311;

/// Codes the choices of every colour with the model.
pub fn encode(
    graph: &StoredGraph,
    choices: &Choices,
) -> Result<Vec<u8>, String> {
    let mut encoder = Encoder::new();
    match graph {
        StoredGraph::Csr(g) => code_all(&mut encoder, g, choices.streams.len(), Some(choices))?,
        StoredGraph::Chains(g) => code_all(&mut encoder, g, choices.streams.len(), Some(choices))?,
    };
    Ok(encoder.finish())
}

/// Decodes the choices of `n_colors` colours from `bytes`, see [encode].
pub fn decode(
    graph: &StoredGraph,
    n_colors: usize,
    bytes: &[u8],
) -> Result<Choices, String> {
    let mut decoder = Decoder::new(bytes);
    let streams = match graph {
        StoredGraph::Csr(g) => code_all(&mut decoder, g, n_colors, None)?,
        StoredGraph::Chains(g) => code_all(&mut decoder, g, n_colors, None)?,
    };
    if decoder.overrun() {
        return Err("choice stream truncated".to_string())
    }
    Ok(Choices { streams })
}

/// A two-layer mixer state for one decision.
struct Mixed {
    x: [i32; N_IN],
    fx: [i32; N_FIN],
    sets: [usize; 3],
    p_sub: [u32; 3],
    p_mix: u32,
    apm: usize,
}

struct Model {
    t: Tables,
    trees: Vec<Vec<u32>>,
    edge_ctr: Vec<u32>,
    copy_ctr: Vec<u32>,
    vote_ctr: Vec<u32>,
    weights: [Vec<i32>; 3],
    final_weights: Vec<i32>,
    apm: Apm,
    /// Visits to every node, and how often every edge was taken.
    visits: Vec<u32>,
    taken: Vec<u32>,
    /// Most recent colours (up to [N_TOOK], newest last) that took every edge.
    took: Vec<[u32; N_TOOK]>,
    /// The most recent colour at every node and its rank.
    last_at: Vec<Option<(u32, u32)>>,
    /// Choices of the finished colours: (node, rank) sorted by node, the last
    /// visit to a node; and of the current colour.
    done: Vec<Vec<(u32, u32)>>,
    current: HashMap<u32, u32>,
}

impl Model {
    fn new<G: GraphWalk>(g: &G) -> Model {
        let t = Tables::new();
        let max_n = (0..g.n_nodes()).map(|j| { let (a, b) = g.row(j); b - a }).max().unwrap_or(1) as u32;
        let max_bits = 32 - max_n.saturating_sub(1).leading_zeros();
        let w0 = (W_ONE * 3 / 10) as i32;
        let apm = Apm::new(&t, N_CLASSES, ApmRate::Fixed(APM_RATE));
        Model {
            t,
            trees: (0..=max_bits).map(|b| vec![C0; 1 << b]).collect(),
            edge_ctr: vec![C0; g.n_edges()],
            copy_ctr: vec![C0; N_KINDS * 8],
            vote_ctr: vec![C0; 8],
            weights: [vec![w0; SETS_A * N_IN], vec![w0; SETS_B * N_IN], vec![w0; SETS_C * N_IN]],
            final_weights: vec![(W_ONE * 33 / 100) as i32; N_CLASSES * N_FIN],
            apm,
            visits: vec![0; g.n_nodes()],
            taken: vec![0; g.n_edges()],
            took: vec![[u32::MAX; N_TOOK]; g.n_edges()],
            last_at: vec![None; g.n_nodes()],
            done: Vec::new(),
            current: HashMap::new(),
        }
    }

    /// The choice of `color` (the current or a finished colour) at `node`.
    fn choice_of(&self, node: u32, color: u32) -> Option<u32> {
        if color as usize == self.done.len() {
            return self.current.get(&node).copied()
        }
        let made = self.done.get(color as usize)?;
        made.binary_search_by_key(&node, |x| x.0).ok().map(|i| made[i].1)
    }

    fn dot(&self, w: &[i32], x: &[i32]) -> i32 {
        let dot: i64 = x.iter().zip(w).map(|(&a, &b)| a as i64 * b as i64).sum();
        ((dot / W_ONE) as i32).clamp(-SMAX, SMAX)
    }

    /// Mixes the inputs `x` with the weight sets `sets` of `class`.
    fn mix(&self, x: [i32; N_IN], sets: [usize; 3], class: usize) -> (Mixed, u32) {
        let mut fx = [BIAS; N_FIN];
        let mut p_sub = [0; 3];
        for m in 0..3 {
            let s = self.dot(&self.weights[m][sets[m] * N_IN..][..N_IN], &x);
            fx[m] = s;
            p_sub[m] = self.t.sq(s);
        }
        let s = self.dot(&self.final_weights[class * N_FIN..][..N_FIN], &fx);
        let p_mix = self.t.sq(s);
        let (pa, apm) = self.apm.refine(class, s);
        let p = ((p_mix as i64 + 3 * pa as i64) / 4).clamp(1, 65535) as u32;
        (Mixed { x, fx, sets, p_sub, p_mix, apm }, p)
    }

    fn learn(&mut self, m: &Mixed, class: usize, bit: u32) {
        let round = 1i64 << 15;
        for k in 0..3 {
            let err = ((bit as i64) << 16) - m.p_sub[k] as i64;
            for (w, &x) in self.weights[k][m.sets[k] * N_IN..][..N_IN].iter_mut().zip(&m.x) {
                *w += ((x as i64 * err * LR + round) >> 16) as i32;
            }
        }
        let err = ((bit as i64) << 16) - m.p_mix as i64;
        for (w, &x) in self.final_weights[class * N_FIN..][..N_FIN].iter_mut().zip(&m.fx) {
            *w += ((x as i64 * err * LR_FINAL + round) >> 16) as i32;
        }
        self.apm.update(m.apm, bit);
    }
}

/// Codes (`choices` Some) or decodes the choices of colours `0..n_colors`.
fn code_all<C: Coder, G: GraphWalk>(
    c: &mut C,
    g: &G,
    n_colors: usize,
    choices: Option<&Choices>,
) -> Result<Vec<ChoiceStream>, String> {
    if n_colors > 0 && g.n_nodes() <= SINK_NODE as usize {
        return Err(format!("{} colours in a graph of {} nodes", n_colors, g.n_nodes()))
    }
    let mut m = Model::new(g);
    let mut streams = Vec::with_capacity(n_colors);
    let mut refs: Vec<(u32, u32)> = Vec::with_capacity(N_REFS);
    let mut cands: Vec<(u32, usize)> = Vec::new();
    let mut votes: Vec<(u32, usize)> = Vec::with_capacity(N_REFS);
    for col in 0..n_colors as u32 {
        let mut reader = match choices {
            Some(choices) => Some(choices.reader(col).ok_or_else(|| format!("no choice stream for colour {}", col))?),
            None => None,
        };
        let mut writer = ChoiceWriter::new();
        refs.clear();
        m.current.clear();
        let mut node = SOURCE_NODE as usize;
        let (mut steps, mut max_steps) = (0_usize, 1_usize);
        loop {
            let (a, b) = g.row(node);
            let n = (b - a) as u32;
            let rank = match n {
                0 => return Err(format!("colour {} has no edge out of node {}", col, node)),
                1 => 0,
                _ => {
                    let actual = match reader.as_mut() {
                        Some(reader) => reader.next(n).ok_or_else(|| format!("choice stream of colour {} exhausted at node {}", col, node))?,
                        None => 0,
                    };
                    let j = node as u32;
                    if refs.is_empty() {
                        if let Some((last, _)) = m.last_at[node] { refs.push((last, SCORE_ONE)); }
                    }
                    // Candidates: the references' predictions, the most recent colour's, the previous colour's.
                    cands.clear();
                    votes.clear();
                    for (i, &(rc, _)) in refs.iter().enumerate() {
                        if let Some(p) = m.choice_of(j, rc) {
                            match votes.iter_mut().find(|v| v.0 == p) {
                                Some(v) => v.1 += 1,
                                None => { votes.push((p, 1)); cands.push((p, i)); },
                            }
                        }
                    }
                    if let Some((_, p)) = m.last_at[node] {
                        if !cands.iter().any(|x| x.0 == p) { cands.push((p, KIND_LAST)); }
                    }
                    if let Some(p) = col.checked_sub(1).and_then(|prev| m.choice_of(j, prev)) {
                        if !cands.iter().any(|x| x.0 == p) { cands.push((p, KIND_PREV)); }
                    }
                    if n <= MAX_LISTED {
                        for r in 0..n {
                            if !cands.iter().any(|x| x.0 == r) { cands.push((r, KIND_RANK)); }
                        }
                    }
                    let nb = (n.min(4) as usize - 2).min(2);
                    let top_score = refs.first().map_or(0, |r| match r.1 / SCORE_ONE { 0 => 0, s => (31 - s.leading_zeros()).min(7) as usize });
                    let vis = (32 - m.visits[node].leading_zeros()).min(11) as usize;
                    let mut result = None;
                    for idx in 0..cands.len() {
                        let (cand, kind) = cands[idx];
                        if n <= MAX_LISTED && idx + 1 == cands.len() {
                            result = Some(cand);
                            break
                        }
                        let edge = a + cand as usize;
                        let is_copy = kind < KIND_RANK;
                        let class = if kind == 0 { 0 } else if kind < N_REFS { 1 } else if is_copy { 2 } else { 3 };
                        let copy_ctx = kind * 8 + idx.min(7);
                        let v = votes.iter().find(|x| x.0 == cand).map_or(0, |x| x.1).min(7);
                        // Share of the visits that took the edge, (taken + 0.4) / (visits + 0.8).
                        let share = ((m.taken[edge] as u64 * 10 + 4) << 16) / (m.visits[node] as u64 * 10 + 8);
                        let x = [
                            m.t.st(m.edge_ctr[edge]),
                            if is_copy { m.t.st(m.copy_ctr[copy_ctx]) } else { 0 },
                            m.t.st(m.vote_ctr[v]),
                            if is_copy { FLAG } else { 0 },
                            m.t.stretch[share.clamp(66, 65470) as usize],
                            BIAS,
                        ];
                        let sets = [class * 3 + nb, (class * 8 + v) * 8 + top_score, (class * 12 + vis) * 4 + idx.min(3)];
                        let (mixed, p) = m.mix(x, sets, class);
                        let bit = c.code((actual == cand) as u32, p);
                        m.learn(&mixed, class, bit);
                        update(&mut m.edge_ctr[edge], bit);
                        if is_copy { update(&mut m.copy_ctr[copy_ctx], bit); }
                        update(&mut m.vote_ctr[v], bit);
                        if bit == 1 {
                            result = Some(cand);
                            break
                        }
                    }
                    let rank = match result {
                        Some(rank) => rank,
                        None => {
                            let bits = 32 - (n - 1).leading_zeros();
                            code_tree(c, &mut m.trees[bits as usize], bits, actual).min(n - 1)
                        },
                    };
                    // Reference scores; on a miss of all, the most recent colour that took the edge joins.
                    let mut any_right = false;
                    for i in 0..refs.len() {
                        if let Some(p) = m.choice_of(j, refs[i].0) {
                            if p == rank { refs[i].1 = refs[i].1.saturating_add(SCORE_ONE); any_right = true; } else { refs[i].1 >>= 1; }
                        }
                    }
                    let edge = a + rank as usize;
                    if !any_right {
                        if let Some(&new) = m.took[edge].iter().rev().find(|&&x| x != u32::MAX && !refs.iter().any(|r| r.0 == x)) {
                            if refs.len() < N_REFS {
                                refs.push((new, SCORE_ONE));
                            } else if let Some(min) = refs.iter_mut().min_by_key(|r| r.1) {
                                *min = (new, SCORE_ONE);
                            }
                        }
                    }
                    refs.sort_by_key(|r| std::cmp::Reverse(r.1));
                    m.current.insert(j, rank);
                    let took = &mut m.took[edge];
                    let keep: Vec<u32> = took.iter().copied().filter(|&x| x != u32::MAX && x != col).collect();
                    let start = keep.len().saturating_sub(N_TOOK - 1);
                    *took = [u32::MAX; N_TOOK];
                    for (slot, &x) in took.iter_mut().zip(keep[start..].iter().chain(std::iter::once(&col))) { *slot = x; }
                    m.last_at[node] = Some((col, rank));
                    m.visits[node] += 1;
                    m.taken[edge] += 1;
                    writer.push(rank, n);
                    rank
                },
            };
            let edge = a + rank as usize;
            if node == SOURCE_NODE as usize {
                max_steps = g.max_steps(edge);
            }
            steps += 1;
            node = g.target(edge);
            if node == SINK_NODE as usize {
                break
            }
            if steps >= max_steps {
                return Err(format!("colour {} did not reach the sink in {} steps", col, max_steps))
            }
            if node >= g.n_nodes() {
                return Err(format!("colour {} reached node {} of {}", col, node, g.n_nodes()))
            }
        }
        if let Some(reader) = reader {
            if !reader.is_finished() {
                return Err(format!("unused choices left for colour {}", col))
            }
        }
        let mut made: Vec<(u32, u32)> = m.current.drain().collect();
        made.sort_unstable();
        m.done.push(made);
        streams.push(writer.finish());
    }
    Ok(streams)
}

#[cfg(test)]
mod tests {
    use std::ops::Range;

    use super::*;
    use super::super::{GraphBuilder, GraphWalk};

    /// Related colours: copies of a few ancestral paths through a sequence of
    /// bubbles, with new variants; starting from one of 40 first entries (a
    /// source node with more than [MAX_LISTED] edges).
    fn dicts() -> Vec<Vec<(usize, Range<usize>)>> {
        let mut x = 0x9E37_79B9_7F4A_7C15_u64;
        let mut rand = move |n: u64| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; (x >> 11) % n };
        let entry = |r: u64| ((1 + r % 5) as usize, (r * 10) as usize..(r * 10 + 1) as usize);
        // Position p is entry p, or one of its variants 1000 + 4p + v.
        let variant = |p: u64, rate: u64, rand: &mut dyn FnMut(u64) -> u64| if rand(rate) == 0 { 1000 + 4 * p + rand(4) } else { p };
        let ancestors: Vec<Vec<u64>> = (0..4).map(|_| (0..80).map(|p| variant(p, 5, &mut rand)).collect()).collect();
        (0..150).map(|i| {
            let mut path = vec![100_000 + i % 40];
            path.extend(ancestors[rand(4) as usize].iter().enumerate().map(|(p, &r)| if rand(30) == 0 { variant(p as u64, 1, &mut rand) } else { r }));
            path.into_iter().map(entry).collect()
        }).collect()
    }

    fn builder(dicts: &[Vec<(usize, Range<usize>)>]) -> GraphBuilder {
        let mut builder = GraphBuilder::new();
        for d in dicts {
            builder.add_path(d.iter().map(|x| x.0).sum(), d);
        }
        builder
    }

    fn roundtrip(graph: &StoredGraph, choices: &Choices) -> usize {
        let bytes = encode(graph, choices).unwrap();
        let decoded = decode(graph, choices.streams.len(), &bytes).unwrap();
        assert_eq!(&decoded, choices);
        bytes.len()
    }

    #[test]
    fn roundtrip_both_layouts() {
        let dicts = dicts();
        let (csr, csr_choices) = builder(&dicts).to_csr();
        let (chains, choices) = builder(&dicts).to_chains();
        let (start, end) = csr.row(SOURCE_NODE as usize);
        assert!(end - start > MAX_LISTED as usize);
        let (csr, chains) = (StoredGraph::Csr(csr), StoredGraph::Chains(chains));

        let n_bits = 8 * choices.streams.iter().map(|s| s.bits.len()).sum::<usize>();
        let coded = roundtrip(&chains, &choices);
        assert_eq!(roundtrip(&csr, &csr_choices), coded);
        // Related colours make the same choices.
        assert!(coded * 8 * 2 < n_bits, "{} bytes for {} bits", coded, n_bits);

        // The decoded paths walk to the same dictionaries.
        let decoded = decode(&chains, dicts.len(), &encode(&chains, &choices).unwrap()).unwrap();
        for color in 0..dicts.len() as u32 {
            assert_eq!(chains.walk_path(&decoded, color).unwrap(), chains.walk_path(&choices, color).unwrap());
        }
    }

    #[test]
    fn no_colours() {
        let (graph, _) = builder(&dicts()).finish();
        assert_eq!(roundtrip(&graph, &Choices::default()), 4);
    }

    #[test]
    fn wrong_choices_and_corrupt_streams() {
        let dicts = dicts();
        let (graph, mut choices) = builder(&dicts).finish();
        let mut bytes = encode(&graph, &choices).unwrap();

        // Fewer colours than coded, more, and a truncated stream.
        assert!(decode(&graph, dicts.len() - 1, &bytes).is_ok());
        assert!(decode(&graph, dicts.len() + 1, &bytes).is_err());
        bytes.truncate(bytes.len() / 2);
        assert!(decode(&graph, dicts.len(), &bytes).is_err());

        // A stream with a choice left over.
        choices.streams[1].bits.extend([0xff; 4]);
        assert!(encode(&graph, &choices).is_err());
        choices.streams.truncate(1);
        assert!(encode(&graph, &choices).is_ok());
    }
}
