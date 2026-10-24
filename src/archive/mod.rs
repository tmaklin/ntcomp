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
//! Compressing collections of assemblies against an index into an NLZ-OG
//! archive ([ArchiveWriter]) and decompressing any of their contigs
//! ([ArchiveReader]). The archive format is [crate::graph]; on top of it the
//! archive keeps a hash, a name and the input file of every contig, and puts
//! back the strand and the edits of the contigs.
//!
pub mod reader;
pub mod writer;

pub use reader::{ArchiveReader, ContigDecoder, DecodedContig};
pub use writer::{ArchiveOpts, ArchiveWriter, InputFile};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::related_index;

    #[test]
    fn archives_roundtrip_with_random_access() {
        let (seqs, index) = related_index(31);
        // Two input files; one contig reverse complemented so the strand is chosen.
        let mut contigs: Vec<(Vec<u8>, Vec<u8>)> = seqs.iter().enumerate().map(|(i, seq)| (format!("contig{} sample", i).into_bytes(), seq.clone())).collect();
        contigs[2].1 = crate::strands::reverse_complement(&contigs[2].1);
        let files: Vec<InputFile> = vec![(b"first".to_vec(), contigs[..1].to_vec()), (b"second".to_vec(), contigs[1..].to_vec())];

        let opts = ArchiveOpts { collapse_chains: true, reorient: true, edit: true, ..Default::default() };
        let mut writer = ArchiveWriter::new(opts, &index);
        assert_eq!(writer.call_edits(&files), vec![0, 1, 2, 3]);
        assert_eq!(writer.add_files(&files).unwrap(), vec![0, 1, 2, 3]);
        let mut bytes = Vec::new();
        writer.write(&mut bytes).unwrap();

        let reader = ArchiveReader::open(&mut bytes.as_slice(), &index).unwrap();
        assert_eq!((reader.n_contigs(), reader.n_files()), (4, 2));
        assert_eq!((reader.file_name(0), reader.file_name(3)), (Some(&b"first"[..]), Some(&b"second"[..])));
        assert_eq!(reader.contig_name(1), Some(&b"contig1 sample"[..]));

        // Every contig, in reverse order, two at a time.
        let colors: Vec<u32> = (0..4).rev().collect();
        let mut decoder = reader.decoder(&colors).unwrap();
        for chunk in colors.chunks(2) {
            for decoded in decoder.decode(chunk).unwrap() {
                let decoded = decoded.unwrap();
                assert_eq!(decoded.sequence, contigs[decoded.color as usize].1, "color {}", decoded.color);
                assert_eq!(decoded.hash_match, Some(true));
            }
        }

        // One contig by name (the whole name or its first word).
        assert_eq!(reader.find_contig("contig2 sample").unwrap(), vec![2]);
        assert_eq!(reader.find_contig("contig2").unwrap(), vec![2]);
        let decoded = reader.decoder(&[2]).unwrap().decode(&[2]).unwrap().pop().unwrap().unwrap();
        assert_eq!(decoded.sequence, contigs[2].1);

        // Unknown names and colours are errors.
        assert!(reader.find_contig("contig9").is_err());
        assert!(reader.decoder(&[4]).is_err());
    }
}
