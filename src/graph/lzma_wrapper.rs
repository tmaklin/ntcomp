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
//! LZMA compression for archive sections.
//!
//! Sections are compressed as raw LZMA streams without the `.lzma` header or
//! an end marker. [LzmaData] stores the parameters and the uncompressed size
//! needed to decompress the stream.
//!
use std::io::{Read, Write};

use lzma_rust2::{LzmaOptions, LzmaReader, LzmaWriter};

type E = Box<dyn std::error::Error>;

/// LZMA preset used for compressing.
const LZMA_PRESET: u32 = 9;

/// Smallest LZMA dictionary used for compressing.
const MIN_DICT_SIZE: u32 = 4096;

#[derive(Debug)]
pub struct LzmaError(pub String);

impl std::fmt::Display for LzmaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "LZMA: {}", self.0)
    }
}

impl std::error::Error for LzmaError {}

/// Raw LZMA stream and the parameters needed to decompress it.
#[derive(serde::Serialize, serde::Deserialize, Debug, Clone, PartialEq, Eq)]
pub struct LzmaData {
    /// Size of the uncompressed data in bytes.
    pub uncompressed_len: u64,

    /// LZMA properties byte (lc, lp and pb).
    pub props: u8,

    /// LZMA dictionary size in bytes.
    pub dict_size: u32,

    /// Raw LZMA stream without a header or an end marker.
    pub bytes: Vec<u8>,
}

/// Compresses `data`.
///
/// Uses no literal context bits (lc = lp = pb = 0): this is best for the
/// bit-packed sections and as good as the defaults for the contig names. The
/// dictionary is not larger than the data.
pub fn compress(
    data: &[u8],
) -> Result<LzmaData, E> {
    let mut options = LzmaOptions::with_preset(LZMA_PRESET);
    options.lc = 0;
    options.lp = 0;
    options.pb = 0;
    let len: u32 = data.len().try_into().unwrap_or(u32::MAX);
    options.dict_size = options.dict_size.min(len.max(MIN_DICT_SIZE));

    let mut writer = LzmaWriter::new_no_header(Vec::new(), &options, false)?;
    writer.write_all(data)?;
    let props = writer.props();
    let bytes = writer.finish()?;

    Ok(LzmaData {
        uncompressed_len: data.len() as u64,
        props,
        dict_size: options.dict_size,
        bytes,
    })
}

/// Decompresses `lzma`; fails if the stream is invalid or does not hold
/// [LzmaData::uncompressed_len] bytes.
pub fn decompress(
    lzma: &LzmaData,
) -> Result<Vec<u8>, E> {
    let mut reader = LzmaReader::new_with_props(lzma.bytes.as_slice(), lzma.uncompressed_len, lzma.props, lzma.dict_size, None)?;
    let mut data = Vec::new();
    reader.read_to_end(&mut data)?;
    if data.len() as u64 != lzma.uncompressed_len {
        return Err(Box::new(LzmaError(format!("expected {} bytes, got {}", lzma.uncompressed_len, data.len()))))
    }
    Ok(data)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lzma_roundtrip() {
        for data in [Vec::new(), b"a".to_vec(), b"OZ038544.1 Escherichia coli isolate 30345_1#313".repeat(100)] {
            let lzma = compress(&data).unwrap();
            assert_eq!(decompress(&lzma).unwrap(), data);
        }
    }

    #[test]
    fn repeated_data_compresses() {
        let data = b"OZ038544.1 Escherichia coli isolate 30345_1#313".repeat(1000);
        assert!(compress(&data).unwrap().bytes.len() * 100 < data.len());
    }

    #[test]
    fn truncated_stream_is_an_error() {
        let data: Vec<u8> = (0..20000_u32).map(|i| (i.wrapping_mul(2654435761) >> 24) as u8).collect();
        let mut lzma = compress(&data).unwrap();
        lzma.bytes.truncate(lzma.bytes.len() / 2);
        assert!(decompress(&lzma).is_err());
    }
}
