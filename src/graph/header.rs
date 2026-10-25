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
//! The NLZ-OG archive header: the format, the index the archive was written
//! with, and the sizes of the sections that follow it.
//!
//! **Layout** ([FILE_HEADER_BYTES] bytes, integers little-endian): the magic
//! bytes [MAGIC], the format [VERSION] (u16), the index's _k_ (u32), number
//! of _k_-mers (u64) and fingerprint ([FINGERPRINT_BYTES] bytes; all three
//! zero for an archive written without an index, see [IndexId]), the number
//! of contigs (u32), then the sizes (u64) of the sections in their order in
//! the file: file ranges, contig names, graph, hashes, strands, SNPs and
//! choices.
//!
use std::io::{Read, Write};

use crate::index::{IndexId, FINGERPRINT_BYTES};

type E = Box<dyn std::error::Error>;

/// The first bytes of an archive.
pub const MAGIC: &[u8; 4] = b"NLOG";

/// The format version this crate writes and reads.
pub const VERSION: u16 = 1;

/// Bytes of the header.
pub const FILE_HEADER_BYTES: usize = 4 + 2 + 4 + 8 + FINGERPRINT_BYTES + 4 + 7 * 8;

/// Header of an archive written by [GraphEncoder](super::encoder::GraphEncoder).
#[derive(Clone, Debug, PartialEq)]
pub struct FileHeader {
    /// The index the archive was written with (None without an index).
    pub index: Option<IndexId>,

    /// Number of contigs (colours).
    pub n_contigs: u32,

    pub file_range_bytes: u64,

    pub contig_name_bytes: u64,

    pub graph_bytes: u64,

    pub hash_bytes: u64,

    /// Size of the section that marks the reverse complemented contigs.
    pub strand_bytes: u64,

    /// Size of the section with the edited SNPs.
    pub snp_bytes: u64,

    pub choice_bytes: u64,
}

#[derive(Debug)]
pub struct HeaderError(pub String);

impl std::fmt::Display for HeaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for HeaderError {}

impl FileHeader {
    /// The section sizes in their order in the file.
    fn section_bytes(&self) -> [u64; 7] {
        [self.file_range_bytes, self.contig_name_bytes, self.graph_bytes, self.hash_bytes, self.strand_bytes, self.snp_bytes, self.choice_bytes]
    }

    /// Total size of the sections, None if it overflows.
    pub fn total_section_bytes(&self) -> Option<u64> {
        self.section_bytes().iter().try_fold(0_u64, |a, &b| a.checked_add(b))
    }

    /// The header's [FILE_HEADER_BYTES] bytes.
    pub fn to_bytes(&self) -> [u8; FILE_HEADER_BYTES] {
        let mut out = Vec::with_capacity(FILE_HEADER_BYTES);
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&VERSION.to_le_bytes());
        let id = self.index.unwrap_or(IndexId { k: 0, n_kmers: 0, fingerprint: [0; FINGERPRINT_BYTES] });
        out.extend_from_slice(&id.k.to_le_bytes());
        out.extend_from_slice(&id.n_kmers.to_le_bytes());
        out.extend_from_slice(&id.fingerprint);
        out.extend_from_slice(&self.n_contigs.to_le_bytes());
        for n in self.section_bytes() {
            out.extend_from_slice(&n.to_le_bytes());
        }
        out.try_into().expect("the header has FILE_HEADER_BYTES bytes")
    }

    /// Writes the header to `out`.
    pub fn write_to<W: Write>(
        &self,
        out: &mut W,
    ) -> std::io::Result<()> {
        out.write_all(&self.to_bytes())
    }

    /// Reads a header from its bytes: an error if they are not an archive
    /// header of [VERSION].
    pub fn from_bytes(
        bytes: &[u8; FILE_HEADER_BYTES],
    ) -> Result<Self, E> {
        if &bytes[..4] != MAGIC {
            if bytes[..6].iter().all(|&b| b == 0) {
                return Err(Box::new(HeaderError("the archive was written by an older version of ntcomp; encode it again".to_string())))
            }
            return Err(Box::new(HeaderError("not an ntcomp archive".to_string())))
        }
        let version = u16::from_le_bytes([bytes[4], bytes[5]]);
        if version != VERSION {
            return Err(Box::new(HeaderError(format!("unsupported archive version {} (this ntcomp reads version {})", version, VERSION))))
        }
        let mut at = 6;
        let mut take = |n: usize| { let field = &bytes[at..at + n]; at += n; field };
        let k = u32::from_le_bytes(take(4).try_into()?);
        let n_kmers = u64::from_le_bytes(take(8).try_into()?);
        let fingerprint: [u8; FINGERPRINT_BYTES] = take(FINGERPRINT_BYTES).try_into()?;
        let index = (k != 0).then_some(IndexId { k, n_kmers, fingerprint });
        let n_contigs = u32::from_le_bytes(take(4).try_into()?);
        let mut sizes = [0_u64; 7];
        for size in sizes.iter_mut() {
            *size = u64::from_le_bytes(take(8).try_into()?);
        }
        let [file_range_bytes, contig_name_bytes, graph_bytes, hash_bytes, strand_bytes, snp_bytes, choice_bytes] = sizes;
        let header = FileHeader { index, n_contigs, file_range_bytes, contig_name_bytes, graph_bytes, hash_bytes, strand_bytes, snp_bytes, choice_bytes };
        if header.total_section_bytes().is_none() {
            return Err(Box::new(HeaderError("corrupt archive header: the section sizes overflow".to_string())))
        }
        Ok(header)
    }

    /// Reads a header from `input`, see [FileHeader::from_bytes].
    pub fn read_from<R: Read>(
        input: &mut R,
    ) -> Result<Self, E> {
        let mut bytes = [0_u8; FILE_HEADER_BYTES];
        input.read_exact(&mut bytes).map_err(|e| HeaderError(format!("cannot read the archive header: {}", e)))?;
        Self::from_bytes(&bytes)
    }

    /// Checks that the archive was written with `index`: an error naming
    /// both if not. An archive written without an index passes.
    pub fn check_index(
        &self,
        index: IndexId,
    ) -> Result<(), E> {
        match self.index {
            Some(stored) if stored != index => Err(Box::new(HeaderError(format!(
                "the archive was written with a different index ({}) than the one given ({})", stored, index)))),
            _ => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> FileHeader {
        FileHeader {
            index: Some(IndexId { k: 31, n_kmers: 1234, fingerprint: [1, 2, 3, 4, 5, 6, 7, 8] }),
            n_contigs: 5,
            file_range_bytes: 1, contig_name_bytes: 2, graph_bytes: 3, hash_bytes: 4, strand_bytes: 5, snp_bytes: 6, choice_bytes: 7,
        }
    }

    #[test]
    fn roundtrip() {
        for h in [header(), FileHeader { index: None, ..header() }] {
            let mut bytes = Vec::new();
            h.write_to(&mut bytes).unwrap();
            assert_eq!(bytes.len(), FILE_HEADER_BYTES);
            assert_eq!(&bytes[..4], MAGIC);
            assert_eq!(FileHeader::read_from(&mut bytes.as_slice()).unwrap(), h);
        }
    }

    #[test]
    fn rejects_other_files_and_versions() {
        let bytes = header().to_bytes();
        let err = |b: [u8; FILE_HEADER_BYTES]| FileHeader::from_bytes(&b).unwrap_err().to_string();
        let mut other = bytes;
        other[0] = b'X';
        assert!(err(other).contains("not an ntcomp archive"));
        let mut old = bytes;
        old[..6].fill(0);
        assert!(err(old).contains("older version"));
        let mut newer = bytes;
        newer[4] = 2;
        assert!(err(newer).contains("version 2"));
        let mut overflow = bytes;
        overflow[FILE_HEADER_BYTES - 16..].fill(0xff);
        assert!(err(overflow).contains("overflow"));
        assert!(FileHeader::read_from(&mut &bytes[..10]).is_err());
    }

    #[test]
    fn index_check() {
        let h = header();
        let id = h.index.unwrap();
        assert!(h.check_index(id).is_ok());
        assert!(h.check_index(IndexId { n_kmers: 1235, ..id }).is_err());
        assert!(h.check_index(IndexId { fingerprint: [0; FINGERPRINT_BYTES], ..id }).is_err());
        assert!(FileHeader { index: None, ..h }.check_index(id).is_ok());
    }
}
