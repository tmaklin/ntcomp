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

use std::collections::{
    HashSet,
};

use std::fs::File;

use blake3::{
    hash,
    Hash,
};

use indexmap::IndexSet;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use petgraph::graph::Graph;

use clap::Parser;
use log::info;
use needletail::Sequence;
use needletail::parser::SequenceRecord;

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
) -> Option<SequenceRecord<'_>> {
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

            // Number of records to store in each block
            let block_size = 65536;

            let header_bytes = ntcomp::encode_file_header(0,0,0,0).unwrap();
            let _ = stdout.write_all(&header_bytes);

            info!("Encoding fastX data...");
            let mut reader = needletail::parse_fastx_file(query_file).unwrap_or_else(|_| panic!("Expected valid fastX file"));
            let mut dictionaries: Vec<Vec<(usize, std::ops::Range<usize>)>> = Vec::new();
            let mut num_records = 0;

            while let Some(rec) = read_from_fastx_parser(&mut *reader) {
                let seqrec = rec.normalize(true);
                num_records += 1;

                dictionaries.push(ntcomp::encode_sequence(&seqrec, &sbwt, &lcs).unwrap());

                if num_records % block_size == 0 {
                    let u64_encoding = dictionaries.iter().flat_map(|x| ntcomp::encode::encode_dictionary(x, &sbwt).unwrap()).collect::<Vec<u64>>();
                    let _ = ntcomp::write_block_to(&u64_encoding, block_size, &mut stdout);
                    dictionaries.clear();
                }
            }
            if num_records % block_size != 0 {
                let u64_encoding = dictionaries.iter().flat_map(|x| ntcomp::encode::encode_dictionary(x, &sbwt).unwrap()).collect::<Vec<u64>>();
                let _ = ntcomp::write_block_to(&u64_encoding, num_records % block_size, &mut stdout);
            }

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
            let _file_header = ntcomp::decode_file_header(&header_bytes).unwrap();

            info!("Decoding encoded data...");
            let mut i = 0;
            while let Ok(encoding) = ntcomp::read_block(&_file_header, &sbwt, &mut conn) {
                let records = ntcomp::decode_sequence(&encoding, &sbwt);
                records.iter().for_each(|nucleotides| {
                    let _ = writeln!(&mut stdout, ">seq.{}", i + 1);
                    let _ = writeln!(&mut stdout,
                                     "{}", nucleotides.iter().map(|x| *x as char).collect::<String>());
                    let _ = stdout.flush();
                    i += 1;
                });
            }
        },
        Some(cli::Commands::View {
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
            let _file_header = ntcomp::decode_file_header(&header_bytes).unwrap();

            info!("Decoding encoded data...");
            while let Ok(encoding) = ntcomp::read_block(&_file_header, &sbwt, &mut conn) {
                let plaintext = ntcomp::decode_positions(&encoding, &sbwt);
                plaintext.iter().for_each(|(contig_id, start, end, length, encoding_type, colex)| {
                    writeln!(&mut stdout, "{}\t{}\t{}\t{}\t{}\t{}", contig_id, start, end, length, encoding_type, colex.unwrap_or(0)).unwrap();
                    stdout.flush().unwrap();
                });
            }
        },
        Some(cli::Commands::Collection {
            query_files,
            index_prefix,
        }) => {
            init_log(2);
            let mut stdout = BufWriter::new(std::io::stdout());

            let (sbwt, lcs) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());

            let header_bytes = ntcomp::encode_file_header(0,0,0,0).unwrap();
            let _ = stdout.write_all(&header_bytes);

            let mut colex_remapping: IndexSet<u32> = IndexSet::new();
            let mut path_starts: Vec<u32> = Vec::with_capacity(query_files.len());
            let mut path_lengths: Vec<u32> = Vec::with_capacity(query_files.len());
            let mut all_lengths: Vec<Vec<u32>> = Vec::with_capacity(query_files.len());
            for (i, query_file) in query_files.iter().enumerate() {
                eprintln!("{}/{}", i + 1, query_files.len());
                let mut reader = needletail::parse_fastx_file(query_file).unwrap_or_else(|_| panic!("Expected valid fastX file"));

                while let Some(rec) = read_from_fastx_parser(&mut *reader) {
                    let seqrec = rec.normalize(true);

                    let dictionary = ntcomp::encode_sequence(&seqrec, &sbwt, &lcs).unwrap();
                    let (path, lengths) = ntcomp::collection::remap_dictionary(dictionary, &mut colex_remapping);
                    let (mut path_start, mut path_length) = ntcomp::collection::get_path_blocks(&path);
                    path_starts.append(&mut path_start);
                    path_lengths.append(&mut path_length);
                    all_lengths.push(lengths);
                }
            }

            ntcomp::collection::write_paths(&path_starts, &path_lengths, &mut stdout).unwrap();
            ntcomp::collection::write_lengths(&all_lengths, &mut stdout).unwrap();
            ntcomp::collection::write_remapping(&colex_remapping, &mut stdout).unwrap();

            let _ = stdout.flush();
        },
        Some(cli::Commands::Graph {
            query_files,
            index_prefix,
            decompress,
            fasta,
            fasta_columns,
        }) => {
            init_log(2);
            let mut stdout = BufWriter::new(std::io::stdout());

            let (sbwt, lcs) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());

            if *decompress {
                assert!(query_files.len() == 1);
                let mut input = File::open(&query_files[0]).unwrap();
                let mut header_bytes: [u8; 42] = [0; 42];
                input.read_exact(&mut header_bytes).unwrap();
                let header = ntcomp::decode_file_header(&header_bytes).unwrap();

                let mut file_range_bytes = vec![0_u8; header.file_range_bytes as usize];
                input.read_exact(&mut file_range_bytes).unwrap();
                let file_ranges: Vec<(Vec<u8>, core::ops::Range<u32>)> = postcard::from_bytes(&file_range_bytes).unwrap();

                let mut contig_name_bytes = vec![0_u8; header.contig_name_bytes as usize];
                input.read_exact(&mut contig_name_bytes).unwrap();
                let contig_names: Vec<Vec<u8>> = postcard::from_bytes(&contig_name_bytes).unwrap();

                let mut graph_bytes = vec![0_u8; header.graph_bytes as usize];
                input.read_exact(&mut graph_bytes).unwrap();
                let csr: ntcomp::graph::Csr<u32> = postcard::from_bytes(&graph_bytes).unwrap();
                let graph = csr.to_petgraph();

                let mut hash_bytes = vec![0_u8; header.hash_bytes as usize];
                input.read_exact(&mut hash_bytes).unwrap();
                let hashes: Vec<Hash> = postcard::from_bytes(&hash_bytes).unwrap();

                let colors: Vec<u32> = (0..header.n_queries).collect();

                let max_visits = header.max_visits;

                if !fasta {
                    eprintln!("accession\tcontig\tdecoded_len\tpath_found");
                }
                let mut file_idx: usize = 0;
                for (idx, seq) in colors.into_iter().enumerate() {
                    if idx as u32 == file_ranges[file_idx].1.end {
                        file_idx += 1;
                    }
                    let nodes = ntcomp::graph::search(
                        &graph,
                        seq,
                        hashes[idx],
                        &sbwt,
                        max_visits,
                    );

                    let file_name = String::from_utf8(file_ranges[file_idx].0.to_vec()).unwrap_or(file_idx.to_string());
                    let contig_name = String::from_utf8(contig_names[idx].to_vec()).unwrap_or(seq.to_string());
                    if let Some(nodes) = nodes {
                        let sequence = ntcomp::graph::decode_sequence(
                            &graph,
                            &nodes,
                            &sbwt,
                        );
                        if !fasta {
                            eprintln!("{}\t{}\t{}\t{}", file_name, contig_name, sequence.len(), true);
                        } else {
                            stdout.write_all(b">").unwrap();
                            stdout.write_all(&contig_names[idx].to_vec()).unwrap();
                            stdout.write_all(b"\n").unwrap();
                            sequence.chunks(*fasta_columns).for_each(|chunk| {
                                stdout.write_all(chunk).unwrap();
                                stdout.write_all(b"\n").unwrap();
                            });
                            stdout.flush().unwrap();
                        }
                    } else {
                        if !fasta {
                            eprintln!("{}\t{}\t{}\t{}", file_name, contig_name, 0, false);
                        }
                    }
                }
            } else {
                let n_queries = query_files.len();

                let progress = ProgressBar::new(n_queries as u64);
                progress.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}").unwrap());

                let mut graph: Graph<u32, ntcomp::graph::ColexGraphEdge> = Graph::new();
                let mut node_indexes: IndexSet<u32> = IndexSet::new();

                // Dummy nodes to track start and end positions
                let start_node = graph.add_node(0_u32);
                let end_node = graph.add_node(0_u32);
                assert!(start_node.index() == 0);
                assert!(end_node.index() == 1);
                node_indexes.insert(0);
                node_indexes.insert(1);

                let mut hashes: Vec<Hash> = Vec::new();

                let mut contig_names: Vec<Vec<u8>> = Vec::with_capacity(n_queries);
                let mut file_ranges: Vec<(Vec<u8>, core::ops::Range<u32>)> = Vec::with_capacity(n_queries);

                let mut color: u32 = 0;
                for query_file in query_files.iter() {
                    let mut reader = needletail::parse_fastx_file(query_file).unwrap_or_else(|_| panic!("Expected valid fastX file"));
                    let start = color;

                    while let Some(rec) = read_from_fastx_parser(&mut *reader) {
                        contig_names.push(rec.id().to_vec());
                        let seqrec = rec.normalize(true);
                        let sequence_length = seqrec.len();

                        let dictionary = ntcomp::encode_sequence(&seqrec, &sbwt, &lcs).unwrap();
                        let n_entries = dictionary.len();
                        ntcomp::graph::insert_edge(
                            &mut graph,
                            &mut node_indexes,
                            color,
                            &(sequence_length, 0..0),
                            &dictionary[0],
                        );

                        for i in 1..n_entries {
                            ntcomp::graph::insert_edge(
                                &mut graph,
                                &mut node_indexes,
                                color,
                                &dictionary[i - 1],
                                &dictionary[i]
                            );
                        }
                        ntcomp::graph::insert_edge(
                            &mut graph,
                            &mut node_indexes,
                            color,
                            &dictionary[dictionary.len() - 1],
                            &(0, 1..1),
                        );

                        hashes.push(hash(&seqrec));
                        color += 1;
                    }
                    let end = color;
                    let filename: Vec<u8> = query_file.file_prefix().unwrap().as_encoded_bytes().to_vec();
                    file_ranges.push((filename, start..end));

                    progress.inc(1_u64);
                }
                progress.finish();
                let max_visits = ntcomp::graph::deduplicate_edges(&mut graph);

                // TODO move this part to an encoding function
                {

                    let csr = ntcomp::graph::Csr::from_petgraph(&graph);
                    let graph_bytes = postcard::to_allocvec(&csr).unwrap();
                    let hash_bytes = postcard::to_allocvec(&hashes).unwrap();

                    let contig_name_bytes = postcard::to_allocvec(&contig_names).unwrap();
                    let file_range_bytes = postcard::to_allocvec(&file_ranges).unwrap();
                    let header = ntcomp::FileHeader{
                        nlz_header: [0_u8; 6],
                        start_node_bytes: 0,
                        n_queries: color,
                        contig_name_bytes: contig_name_bytes.len().try_into().unwrap(),
                        graph_bytes: graph_bytes.len() as u64,
                        hash_bytes: hash_bytes.len().try_into().unwrap(),
                        file_range_bytes: file_range_bytes.len().try_into().unwrap(),
                        max_visits,
                    };
                    let nbytes = bincode::encode_into_std_write(
                        &header,
                        &mut stdout,
                        bincode::config::standard().with_fixed_int_encoding(),
                    ).unwrap();
                    assert_eq!(nbytes, 42);

                    stdout.write_all(&file_range_bytes).unwrap();
                    stdout.write_all(&contig_name_bytes).unwrap();
                    stdout.write_all(&graph_bytes).unwrap();
                    stdout.write_all(&hash_bytes).unwrap();
                }

                // ntcomp::graph::encode_to(
                //     &graph,
                //     &colex_remapping,
                //     &mut stdout,
                // ).unwrap();
            }

            let _ = stdout.flush();
        },
        None => {},
    }
}
