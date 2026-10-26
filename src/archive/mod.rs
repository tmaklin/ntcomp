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

pub use reader::{ArchiveReader, ContigDecoder, DecodedContig, FileCheck};
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

        let opts = ArchiveOpts::default();
        let mut writer = ArchiveWriter::new(opts, &index);
        assert_eq!(writer.add_files(&files).unwrap(), vec![0, 1, 2, 3]);
        let mut bytes = Vec::new();
        writer.write(&mut bytes).unwrap();

        let reader = ArchiveReader::open(&mut bytes.as_slice(), &index).unwrap();
        assert_eq!((reader.n_contigs(), reader.n_files()), (4, 2));
        assert_eq!((reader.file_name(0), reader.file_name(3)), (Some(&b"first"[..]), Some(&b"second"[..])));
        assert_eq!(reader.contig_name(1), Some(&b"contig1 sample"[..]));

        // Every contig, in reverse order, two at a time: the second file
        // (colours 1 to 3) is checked when its last contig is decoded, over
        // two calls, then the first.
        let colors: Vec<u32> = (0..4).rev().collect();
        let mut decoder = reader.decoder(&colors).unwrap();
        let mut checks = Vec::new();
        for chunk in colors.chunks(2) {
            for decoded in decoder.decode(chunk).unwrap() {
                let decoded = decoded.unwrap();
                assert_eq!(decoded.sequence, contigs[decoded.color as usize].1, "color {}", decoded.color);
            }
            checks.push(decoder.finished_files());
        }
        assert_eq!(checks, vec![vec![], vec![FileCheck { file: 1, hash_match: true }, FileCheck { file: 0, hash_match: true }]]);

        // One contig by name (the whole name or its first word): its file is
        // not checked. Its whole file by the file's name is.
        assert_eq!(reader.find_contig("contig2 sample").unwrap(), vec![2]);
        assert_eq!(reader.find_contig("contig2").unwrap(), vec![2]);
        let mut decoder = reader.decoder(&[2]).unwrap();
        let decoded = decoder.decode(&[2]).unwrap().pop().unwrap().unwrap();
        assert_eq!(decoded.sequence, contigs[2].1);
        assert!(decoder.finished_files().is_empty());
        let colors = reader.find_file("second").unwrap();
        assert_eq!(colors, vec![1, 2, 3]);
        let mut decoder = reader.decoder(&colors).unwrap();
        decoder.decode(&colors).unwrap();
        assert_eq!(decoder.finished_files(), vec![FileCheck { file: 1, hash_match: true }]);

        // A changed contig name changes its file's hash only; a wrong hash
        // fails its file only.
        let mut renamed = files.clone();
        renamed[1].1[0].0 = b"another name".to_vec();
        let mut writer = ArchiveWriter::new(ArchiveOpts::default(), &index);
        writer.add_files(&renamed).unwrap();
        let mut other = Vec::new();
        writer.write(&mut other).unwrap();
        let mut other = ArchiveReader::open(&mut other.as_slice(), &index).unwrap();
        let (hashes, other_hashes) = (reader.hashes.clone().unwrap(), other.hashes.clone().unwrap());
        assert_eq!(hashes[0], other_hashes[0]);
        assert_ne!(hashes[1], other_hashes[1]);
        other.hashes = Some(hashes);
        let all: Vec<u32> = (0..4).collect();
        let mut decoder = other.decoder(&all).unwrap();
        decoder.decode(&all).unwrap();
        assert_eq!(decoder.finished_files(), vec![FileCheck { file: 0, hash_match: true }, FileCheck { file: 1, hash_match: false }]);

        // Unknown names and colours are errors.
        assert!(reader.find_contig("contig9").is_err());
        assert!(reader.find_file("third").is_err());
        assert!(reader.decoder(&[4]).is_err());
    }

    #[test]
    fn runs_of_n_are_cut_out_and_put_back() {
        let (seqs, index) = related_index(31);
        // Runs at the start, in the middle (one long, one a single N, two
        // close together), at the end, a contig of only N, and an IUPAC code.
        let with_runs = |seq: &[u8], runs: &[(usize, usize)]| {
            let mut out = seq.to_vec();
            for &(pos, len) in runs.iter().rev() {
                out.splice(pos..pos, std::iter::repeat_n(b'N', len));
            }
            out
        };
        let contigs: Vec<(Vec<u8>, Vec<u8>)> = vec![
            (b"a".to_vec(), with_runs(&seqs[0], &[(0, 5), (700, 250), (1500, 1), (1510, 3)])),
            (b"b".to_vec(), with_runs(&seqs[1], &[(900, 250), (seqs[1].len(), 10)])),
            (b"c".to_vec(), seqs[2].clone()),
            (b"d".to_vec(), b"NNNN".to_vec()),
            (b"e".to_vec(), with_runs(&seqs[2], &[(300, 1)]).into_iter().enumerate().map(|(i, b)| if i == 300 { b'R' } else { b }).collect()),
        ];
        // Names with the contig's length (whole, with its runs).
        let contigs: Vec<(Vec<u8>, Vec<u8>)> = contigs.into_iter().enumerate()
            .map(|(i, (_, seq))| (format!("NODE_{}_length_{}_cov_1.5", i + 1, seq.len()).into_bytes(), seq)).collect();
        let files: Vec<InputFile> = vec![(b"file".to_vec(), contigs.clone())];
        let mut writer = ArchiveWriter::new(ArchiveOpts::default(), &index);
        writer.add_files(&files).unwrap();
        let mut bytes = Vec::new();
        let summary = writer.write(&mut bytes).unwrap();
        assert!(summary.contains("(8 runs)"), "{}", summary);

        // The lengths the names are coded by are the contigs' lengths.
        let sections = crate::graph::decoder::GraphDecoder::new(&mut bytes.as_slice()).unwrap().decode_from_read_with_index(&index).unwrap();
        let lengths = crate::graph::names::contig_lengths(&sections.graph, &sections.choices, sections.snp_edits.as_ref(), sections.n_runs.as_ref()).unwrap();
        assert_eq!(lengths, contigs.iter().map(|(_, s)| s.len() as u64).collect::<Vec<_>>());

        let reader = ArchiveReader::open(&mut bytes.as_slice(), &index).unwrap();
        let colors: Vec<u32> = (0..5).collect();
        let mut decoder = reader.decoder(&colors).unwrap();
        for decoded in decoder.decode(&colors).unwrap() {
            let decoded = decoded.unwrap();
            assert_eq!(decoded.sequence, contigs[decoded.color as usize].1, "colour {}", decoded.color);
            assert_eq!(reader.contig_name(decoded.color), Some(contigs[decoded.color as usize].0.as_slice()));
        }

        assert_eq!(decoder.finished_files(), vec![FileCheck { file: 0, hash_match: true }]);

        // Without the cut, a parse fails instead of looping.
        assert!(crate::nlz::parse(&contigs[0].1, &index).is_err());
    }
}
