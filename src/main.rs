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
use std::io::{Read, Write};
use std::ops::Range;
use std::path::PathBuf;

use bitpacking::{BitPacker8x, BitPacker};
use bincode::encode_into_std_write;
use bincode::decode_from_std_read;
use bincode::{Encode, Decode};
use bincode::decode_from_slice;
use clap::Parser;
use log::info;
use needletail::Sequence;
use needletail::parser::SequenceRecord;
use sbwt::SbwtIndexVariant;
use sbwt::StreamingIndex;

use dsi_bitstream::prelude::*;

// Command-line interface
mod cli;

#[derive(Encode, Decode)]
struct BlockHeader {
    pub block_size: u32,
    pub num_records: u32,
    pub num_u64: u32,
    pub encoded_size: u32,
    pub rice_param: u8,
    pub bitpacker_exponent: u8,
    pub placeholder1: u64,
    pub placeholder2: u32,
    pub placeholder3: u16,
    // TODO should store the used codec so others can be used
    // TODO store total number of bases?
}

/// Initializes the logger with verbosity given in `log_max_level`.
fn init_log(log_max_level: usize) {
    stderrlog::new()
    .module(module_path!())
    .quiet(false)
    .verbosity(log_max_level)
    .timestamp(stderrlog::Timestamp::Off)
    .init()
    .unwrap();
}

// Given a needletail parser, reads the next contig sequence
fn read_from_fastx_parser(
    reader: &mut dyn needletail::parser::FastxReader,
) -> Option<SequenceRecord> {
    let rec = reader.next();
    match rec {
        Some(Ok(seqrec)) => {
            Some(seqrec)
        },
        None => None,
        Some(Err(_)) => todo!(),
    }
}

// Reads all sequence data from a fastX file
fn read_fastx_file(
    file: &str,
) -> Vec<(String, Vec<u8>)> {
    let mut seq_data: Vec<(String, Vec<u8>)> = Vec::new();
    let mut reader = needletail::parse_fastx_file(file).unwrap_or_else(|_| panic!("Expected valid fastX file at {}", file));
    while let Some(rec) = read_from_fastx_parser(&mut *reader) {
        let seqrec = rec.normalize(true);
        let seqname = String::from_utf8(rec.id().to_vec()).expect("UTF-8");
        seq_data.push((seqname, seqrec.to_vec()));
    }
    seq_data
}

fn read_input_list(
    input_list_file: &String,
    delimiter: u8
) -> Vec<(String, PathBuf)> {
    let fs = match std::fs::File::open(input_list_file) {
        Ok(fs) => fs,
        Err(e) => panic!("  Error in reading --input-list: {}", e),
    };

    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(false)
        .from_reader(fs);

    reader.records().map(|line| {
        if let Ok(record) = line {
            if record.len() > 1 {
                (record[0].to_string(), PathBuf::from(record[1].to_string()))
            } else {
                (record[0].to_string(), PathBuf::from(record[0].to_string()))
            }
        } else {
            panic!("  Error in reading --input-list: {}", input_list_file);
        }
    }).collect::<Vec<(String, PathBuf)>>()
}

fn main() {
    let cli = cli::Cli::parse();

    // Subcommands:
    match &cli.command {
        // copy-paste from kbo-cli
        Some(cli::Commands::Build {
            seq_files,
            input_list,
            output_prefix,
            kmer_size,
            prefix_precalc,
            dedup_batches,
            num_threads,
            mem_gb,
            temp_dir,
            verbose,
        }) => {
            init_log(if *verbose { 2 } else { 1 });

            let mut sbwt_build_options = kbo::BuildOpts::default();
            sbwt_build_options.k = *kmer_size;
            sbwt_build_options.num_threads = *num_threads;
            sbwt_build_options.prefix_precalc = *prefix_precalc;
            sbwt_build_options.dedup_batches = *dedup_batches;
            sbwt_build_options.mem_gb = *mem_gb;
            sbwt_build_options.temp_dir = temp_dir.clone();
            sbwt_build_options.add_revcomp = true;
            sbwt_build_options.build_select = true;

            let mut in_files = seq_files.clone();
            if let Some(list) = input_list {
                let contents = read_input_list(list, b'\t');
                let contents_iter = contents.iter().map(|(_, path)| path.to_str().unwrap().to_string());
                in_files.extend(contents_iter);
            }

            info!("Building SBWT index from {} files...", in_files.len());
            let mut seq_data: Vec<Vec<u8>> = Vec::new();
            in_files.iter().for_each(|file| {
                seq_data.append(&mut read_fastx_file(file).into_iter().map(|(_, seq)| seq).collect::<Vec<Vec<u8>>>());
            });

            let (sbwt, lcs) = kbo::build(&seq_data, sbwt_build_options);

            info!("Serializing SBWT index to {}.sbwt ...", output_prefix.as_ref().unwrap());
            info!("Serializing LCS array to {}.lcs ...", output_prefix.as_ref().unwrap());
            kbo::index::serialize_sbwt(output_prefix.as_ref().unwrap(), &sbwt, &lcs);

        },
        Some(cli::Commands::Encode {
            query_file,
            index_prefix,
        }) => {
            init_log(2);
            let stdout = std::io::stdout();
            info!("Loading SBWT index...");

            let (sbwt, lcs) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());

            // TODO should use MB here
            let block_size = 65536;

            info!("Encoding fastX data...");
            match sbwt {
                SbwtIndexVariant::SubsetMatrix(sbwt) => {
                    let index = StreamingIndex::new(&sbwt, &lcs);

                    let mut reader = needletail::parse_fastx_file(query_file).unwrap_or_else(|_| panic!("Expected valid fastX file"));
                    let mut u64_encoding: Vec<u64> = Vec::new();
                    let mut num_records = 0;

                    while let Some(rec) = read_from_fastx_parser(&mut *reader) {
                        let seqrec = rec.normalize(true);
                        let res: Vec<(usize, Range<usize>)> = index.matching_statistics(&seqrec);
                        num_records += 1;

                        let mut bases = 0;
                        res.iter().rev().for_each(|matches| {
                            if bases == 0 {
                                bases = matches.0;
                                let mut arr: [u8; 8] = [0; 8];
                                arr[0..4].copy_from_slice(&(matches.1.start as u32).to_ne_bytes());
                                arr[4..8].copy_from_slice(&(matches.0 as u32).to_ne_bytes());
                                u64_encoding.push(u64::from_ne_bytes(arr));
                            }
                            bases -= 1;
                        });

                        if u64_encoding.len() > block_size {

                            // Best: Rice
                            let inv_mean: f64 = ((u64_encoding.len() as f64).ln() - (u64_encoding.iter().sum::<u64>() as f64).ln()).exp();
                            let param = dsi_bitstream::codes::rice::log2_b(inv_mean);

                            let word_write = MemWordWriterVec::new(Vec::<u64>::new());
                            let mut writer = BufBitWriter::<BE, _>::new(word_write);

                            u64_encoding.iter().for_each(|n| { writer.write_rice(*n, param).unwrap(); } );
                            let compressed_data = writer.into_inner().unwrap().into_inner();

                            let mut bits = compressed_data.iter().flat_map(|x| {
                                let arr: [u8; 8] = x.to_ne_bytes();
                                let mut arr1: [u8; 4] = [0; 4];
                                let mut arr2: [u8; 4] = [0; 4];
                                arr1[0..4].copy_from_slice(&(arr[0..4]));
                                arr2[0..4].copy_from_slice(&(arr[4..8]));
                                [u32::from_ne_bytes(arr1), u32::from_ne_bytes(arr2)]
                            }).collect::<Vec<u32>>();

                            let mut bitpacked: Vec<u8> = Vec::new();

                            let bitpacker = BitPacker8x::new();
                            let block_size = 256;

                            if block_size - bits.len() % block_size != 0 {
                                bits.resize(bits.len() + (block_size - bits.len() % block_size), 0);
                            }
                            bits.chunks(block_size).for_each(|chunk| {
                                let num_bits: u8 = bitpacker.num_bits(chunk);
                                let mut compressed = vec![0u8; 4 * BitPacker8x::BLOCK_LEN];
                                let clen = bitpacker.compress(chunk, &mut compressed[..], num_bits);

                                assert_eq!(clen, 1024);
                                assert_eq!(num_bits, 32);

                                bitpacked.append(&mut compressed);
                            });

                            let block_header = BlockHeader{ block_size: bitpacked.len() as u32,
                                                            num_records: num_records as u32,
                                                            num_u64: u64_encoding.len() as u32,
                                                            encoded_size: compressed_data.len() as u32,
                                                            rice_param: param as u8,
                                                            bitpacker_exponent: 8_u8,
                                                            placeholder1: 0, placeholder2: 0, placeholder3: 0,
                            };

                            let nbytes = encode_into_std_write(
                                &block_header,
                                &mut stdout.lock(),
                                bincode::config::standard().with_fixed_int_encoding(),
                            );
                            assert_eq!(nbytes.unwrap(), 32);
                            let _ = stdout.lock().write_all(&bitpacked);

                            u64_encoding.clear();
                        }
                        // TODO Write whatever is left in the block
                    }

                },
            };
        },
        Some(cli::Commands::Decode {
            input_path,
            index_prefix,
        }) => {
            init_log(2);
            let stdout = std::io::stdout();
            info!("Loading SBWT index...");
            let (sbwt, _) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());

            info!("Reading encoded data...");
            let mut conn = std::fs::File::open(input_path).unwrap();
            let in_data: Vec<Vec<(u8, usize)>> = bincode::decode_from_std_read(&mut conn, bincode::config::standard()).unwrap();

            info!("Decoding encoded data...");

            match sbwt {
                SbwtIndexVariant::SubsetMatrix(sbwt) => {
                    let k = sbwt.k();

                    let mut i = 1;
                    in_data.iter().for_each(|encoding| {
                        let mut nucleotides: Vec<u8> = Vec::new();

                        encoding.iter().for_each(|(suffix_len, colex)| {
                            let kmer = sbwt.access_kmer(*colex);
                            nucleotides.extend(kmer[(k - (*suffix_len as usize))..k].iter().rev());
                        });

                        let _ = writeln!(&mut stdout.lock(), ">seq.{}", i);
                        let _ = writeln!(&mut stdout.lock(),
                                         "{}", nucleotides.iter().rev().map(|x| *x as char).collect::<String>());

                        i += 1;
                    });
                },
            }
        },
        None => {},
    }
}
