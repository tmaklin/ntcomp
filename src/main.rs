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
use std::io::BufWriter;
use std::io::{Read, Write};
use std::path::PathBuf;

use clap::Parser;
use log::info;
use needletail::Sequence;
use needletail::parser::SequenceRecord;

use ntcomp::decode;
use ntcomp::encode;
use ntcomp::encode_file_header;
use ntcomp::decode_file_header;
use ntcomp::decode_block_header;
use ntcomp::decode_sequence;
use ntcomp::encode_sequence;

mod cli;

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
            let mut stdout = BufWriter::new(std::io::stdout());
            info!("Loading SBWT index...");

            let (sbwt, lcs) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());

            // TODO should use MB here
            let block_size = 65536;

            let header_bytes = encode_file_header(0,0,0,0);
            let _ = stdout.write_all(&header_bytes);

            info!("Encoding fastX data...");
            let mut reader = needletail::parse_fastx_file(query_file).unwrap_or_else(|_| panic!("Expected valid fastX file"));
            let mut u64_encoding: Vec<u64> = Vec::new();
            let mut num_records = 0;

            while let Some(rec) = read_from_fastx_parser(&mut *reader) {
                let seqrec = rec.normalize(true);
                num_records += 1;
                let mut encoding = encode_sequence(&seqrec, &sbwt, &lcs);
                u64_encoding.append(&mut encoding);

                if u64_encoding.len() > block_size {
                    let block = encode::compress_block(&u64_encoding, num_records);
                    let _ = stdout.write_all(&block);
                    num_records = 0;
                    u64_encoding.clear();
                }
            }
            let block = encode::compress_block(&u64_encoding, num_records);
            let _ = stdout.write_all(&block);
            let _ = stdout.flush();
            // TODO rewind back to start and fill file header
        },

        Some(cli::Commands::Decode {
            input_path,
            index_prefix,
        }) => {
            init_log(2);
            let mut stdout = BufWriter::new(std::io::stdout());
            // info!("Loading SBWT index...");
            let (sbwt, _) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());

            // info!("Reading encoded data...");
            let mut conn = std::fs::File::open(input_path).unwrap();

            // File header
            let mut header_bytes: [u8; 32] = [0_u8; 32];
            let _ = conn.read_exact(&mut header_bytes);
            let _file_header = decode_file_header(&header_bytes);

            let mut i = 1;
            while conn.read_exact(&mut header_bytes).is_ok() {
                let header = decode_block_header(&header_bytes);
                let mut bytes: Vec<u8> = vec![0; header.block_size as usize];

                let _ = conn.read_exact(&mut bytes);

                let decompressed = decode::decompress_block(&bytes, &header);
                let decoded = decode::decode_dictionary(&decompressed);

                info!("Decoding encoded data...");
                let mut pointer = 0;
                while pointer < decoded.len() {
                    let (nucleotides, new_pointer) = decode_sequence(&decoded, &sbwt, pointer);
                    pointer = new_pointer;
                    let _ = writeln!(&mut stdout, ">seq.{}", i);
                    let _ = writeln!(&mut stdout,
                                     "{}", nucleotides.iter().rev().map(|x| *x as char).collect::<String>());
                    i += 1;
                }
                let _ = stdout.flush();
            }
        },
        None => {},
    }
}
