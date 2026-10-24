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
//! Reading and writing FASTA and FASTQ files.
//!
use std::io::Write;
use std::path::Path;

use needletail::Sequence;

type E = Box<dyn std::error::Error>;

/// A record: its name and its normalized sequence.
pub type Contig = (Vec<u8>, Vec<u8>);

/// The records of a FASTA or FASTQ file (optionally compressed): each
/// record's name (the header line after `>` or `@`) and its normalized
/// sequence (upper case, other characters than A, C, G and T as N).
pub fn read_contigs<P: AsRef<Path>>(
    path: P,
) -> Result<Vec<Contig>, E> {
    let path = path.as_ref();
    let error = |e: needletail::errors::ParseError| format!("{}: {}", path.display(), e);
    let mut reader = needletail::parse_fastx_file(path).map_err(error)?;
    let mut contigs = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(error)?;
        contigs.push((record.id().to_vec(), record.normalize(true).to_vec()));
    }
    Ok(contigs)
}

/// The sequences of an input to an index: tigs written by `ntcomp tigs`
/// (`.tigs.tac`, see [crate::tigs::read_tigs]) or a FASTA or FASTQ file.
pub fn read_sequences<P: AsRef<Path>>(
    path: P,
) -> Result<Vec<Vec<u8>>, E> {
    let path = path.as_ref();
    if path.to_string_lossy().ends_with(".tigs.tac") {
        let reader = std::io::BufReader::new(std::fs::File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?);
        return crate::tigs::read_tigs(reader).map_err(|e| format!("{}: {}", path.display(), e).into())
    }
    Ok(read_contigs(path)?.into_iter().map(|(_, sequence)| sequence).collect())
}

/// Writes a FASTA record of `name` and `sequence`, `columns` bases per
/// line (at least one).
pub fn write_fasta<W: Write>(
    out: &mut W,
    name: &[u8],
    sequence: &[u8],
    columns: usize,
) -> std::io::Result<()> {
    out.write_all(b">")?;
    out.write_all(name)?;
    out.write_all(b"\n")?;
    for line in sequence.chunks(columns.max(1)) {
        out.write_all(line)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}
