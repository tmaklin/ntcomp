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
use core::ops::Range;

use std::io::Write;

use indexmap::IndexSet;

type E = Box<dyn std::error::Error>;

pub fn remap_dictionary(
    dictionary: Vec<(usize, Range<usize>)>,
    colex_remapping: &mut IndexSet<u32>,
) -> (Vec<u32>, Vec<u32>) {
    let mut path: Vec<u32> = Vec::new();
    let mut lengths: Vec<u32> = Vec::new();
    dictionary.into_iter().for_each(|(length, colex_range)| {
        let colex_rank: u32 = colex_range.start.try_into().unwrap();

        colex_remapping.insert(colex_rank);
        let node_idx: u32 = colex_remapping.get_index_of(&colex_rank).unwrap().try_into().unwrap();
        path.push(node_idx);
        lengths.push(length as u32);
    });
    (path, lengths)
}

pub fn get_path_blocks(
    path: &[u32],
) -> (Vec<u32>, Vec<u32>) {
    let mut prev: u32 = path[0];
    let mut starts: Vec<u32> = vec![prev];
    let mut lengths: Vec<u32> = vec![1];
    let mut i: usize = 0;
    path.iter().skip(1).for_each(|node_idx| {
        if *node_idx == prev + 1 {
            lengths[i] += 1;
        } else {
            starts.push(*node_idx);
            lengths.push(1);
            i += 1;
        }
        prev = *node_idx;
    });
    (starts, lengths)
}

pub fn write_paths<W: Write>(
    starts: &[u32],
    lengths: &[u32],
    out: &mut W,
) -> Result<(), E> {
    let _ = bincode::encode_into_std_write(
        starts,
        out,
        bincode::config::standard(),
    )?;

    let _ = bincode::encode_into_std_write(
        lengths,
        out,
        bincode::config::standard(),
    )?;

    out.flush()?;
    Ok(())
}

pub fn write_lengths<W: Write>(
    all_lengths: &[Vec<u32>],
    out: &mut W,
) -> Result<(), E> {
    let _ = bincode::encode_into_std_write(
        all_lengths,
        out,
        bincode::config::standard(),
    )?;
    out.flush()?;
    Ok(())
}

pub fn write_remapping<W: Write>(
    colex_remapping: &IndexSet<u32>,
    out: &mut W,
) -> Result<(), E> {
    let colex_ranks: Vec<u64> = colex_remapping.iter().map(|colex_rank| *colex_rank as u64).collect::<Vec<u64>>();
    let colex_bytes = crate::encode::minimal_binary_encode(&colex_ranks)?.0;

    // This is already packed so no need to pass through bincode
    let remapping_bytes = colex_bytes.iter().flat_map(|x| {
        x.to_le_bytes()
    }).collect::<Vec<u8>>();
    out.write_all(&remapping_bytes)?;
    out.flush()?;
    Ok(())
}
