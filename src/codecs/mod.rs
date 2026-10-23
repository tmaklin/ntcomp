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
//! Codecs shared by the archive format ([crate::graph]) and the index tigs
//! ([crate::tigs]). None of them depends on the graph or the index.
//!
//! - [arith]: a binary arithmetic coder with adaptive bit models and helpers
//!   for coding numbers with them; used by the archive format and the tig
//!   coder.
//! - [mix]: the stretch and squash tables of logistic mixing, computed in
//!   integer arithmetic; used by the context-mixing models.
//! - [bits]: plain bit-level codes (fixed width, Elias gamma, truncated
//!   binary); used for the colex ranks and the choice streams.
//! - [lzma]: raw LZMA streams for archive sections (contig names, file
//!   ranges).
//! - [tigs]: the context-mixing DNA coder that stores the index tigs
//!   ([crate::tigs::write_tigs], [crate::tigs::read_tigs]).
pub mod arith;
pub mod bits;
pub mod lzma;
pub mod mix;
pub mod tigs;
