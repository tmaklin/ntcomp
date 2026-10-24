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
use ntcomp::graph::GraphWalk;

use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;

use std::fs::File;

use ntcomp::graph::{contig_hash, ContigHash};

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use clap::Parser;
use log::info;
use needletail::Sequence;
use needletail::parser::SequenceRecord;

use rayon::prelude::*;

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

/// Contig names and normalized sequences of `query_file`.
fn read_contigs(
    query_file: &PathBuf,
) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut reader = needletail::parse_fastx_file(query_file).unwrap_or_else(|_| panic!("Expected valid fastX file"));
    let mut contigs = Vec::new();
    while let Some(rec) = read_from_fastx_parser(&mut *reader) {
        contigs.push((rec.id().to_vec(), rec.normalize(true).to_vec()));
    }
    contigs
}

/// Decodes the stored (edited) sequence of `color` from the graph, before the
/// SNPs are put back.
fn decode_edited(
    color: u32,
    graph: &ntcomp::graph::StoredGraph,
    choices: &ntcomp::graph::Choices,
    sbwt: &sbwt::sbwt_index_variant::SbwtIndexVariant,
    lcs: &sbwt::LcsArray,
) -> Result<Vec<u8>, String> {
    let sbwt::sbwt_index_variant::SbwtIndexVariant::SubsetMatrix(index) = sbwt;

    let (sequence_length, dictionary) = graph.walk_path(choices, color).map_err(|e| e.to_string())?;
    let k = index.k() as u32;
    let dictionary = dictionary.into_iter().map(|(colex_rank, len)| (colex_rank, ntcomp::nlz::SuffixLen::from_stored(len, k))).collect();

    ntcomp::nlz::decode_dictionary(dictionary, sequence_length as usize, sbwt, Some(lcs)).map_err(|e| e.to_string())
}

/// Contexts of SNPs in edited sequences, by colour and position, for putting
/// back repeated SNPs (see [ntcomp::graph::snps]).
type SnpContexts = std::collections::HashMap<(u32, u32), [u8; ntcomp::graph::snps::CONTEXT]>;

/// Puts the SNPs of `color` back into its edited sequence, restores its
/// strand and checks its hash. Returns the sequence and whether the hash
/// matches (None if the archive has no hashes).
fn finish_color(
    color: u32,
    mut sequence: Vec<u8>,
    hashes: &Option<Vec<ContigHash>>,
    reversed: &Option<Vec<bool>>,
    snp_edits: &Option<ntcomp::graph::snps::SnpEdits>,
    contexts: &SnpContexts,
) -> Result<(Vec<u8>, Option<bool>), String> {
    let idx = color as usize;
    // Edited SNPs are put back on the stored strand.
    if let Some(snp_edits) = snp_edits {
        snp_edits.revert(color, &mut sequence, &|f, pos| contexts.get(&(f, pos)).copied()).map_err(|e| e.to_string())?;
    }
    // Contigs stored on the other strand are reverse complemented back.
    let sequence = if reversed.as_ref().is_some_and(|r| r[idx]) {
        ntcomp::strands::reverse_complement(&sequence)
    } else {
        sequence
    };

    // Not checked if the archive has no hashes.
    let hash_match = hashes.as_ref().map(|hashes| contig_hash(&sequence) == hashes[idx]);
    Ok((sequence, hash_match))
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
                if file.ends_with(".tigs.tac") {
                    // index tigs written by `ntcomp tigs`
                    let reader = std::io::BufReader::new(File::open(file).unwrap_or_else(|_| panic!("Expected tigs at {}", file)));
                    seq_data.append(&mut ntcomp::tigs::read_tigs(reader).unwrap_or_else(|e| panic!("Could not read tigs from {}: {}", file, e)));
                } else {
                    seq_data.append(&mut read_fastx_file(file).into_iter().map(|(_, seq)| seq).collect::<Vec<Vec<u8>>>());
                }
            });

            let (sbwt, lcs) = kbo::build(&seq_data, sbwt_build_options);

            info!("Serializing SBWT index to {}.sbwt ...", output_prefix.as_ref().unwrap());
            info!("Serializing LCS array to {}.lcs ...", output_prefix.as_ref().unwrap());
            kbo::index::serialize_sbwt(output_prefix.as_ref().unwrap(), &sbwt, &lcs);

        },
        Some(cli::Commands::Tigs {
            seq_files,
            input_list,
            output_prefix,
            fasta,
            index,
            kmer_size,
            min_assemblies,
            gap_join,
            gap_one_sided,
            gap_bubble,
            gap_tip,
            gap_keep_repairable,
            orient_sampling,
            partitions,
            prefix_precalc,
            mem_gb,
            temp_dir,
            num_threads,
            verbose,
        }) => {
            init_log(if *verbose { 2 } else { 1 });

            let mut in_files = seq_files.clone();
            if let Some(list) = input_list {
                let contents = read_input_list(list, b'\t');
                in_files.extend(contents.iter().map(|(_, path)| path.to_str().unwrap().to_string()));
            }

            // The indexes are built on disk: in memory, sorting the k-mers takes
            // k / 4 bytes each.
            let (temp_dir, own_temp_dir) = match temp_dir {
                Some(dir) => (dir.clone(), false),
                None => (format!("{}.tmp", output_prefix), true),
            };
            std::fs::create_dir_all(&temp_dir).unwrap_or_else(|e| panic!("Could not create {}: {}", temp_dir, e));

            let mut opts = ntcomp::tigs::TigOpts::default();
            opts.min_assemblies = *min_assemblies;
            opts.partitions = *partitions;
            opts.orient_sampling = *orient_sampling;
            opts.gaps = ntcomp::tigs::GapRules {
                join: *gap_join,
                one_sided: *gap_one_sided,
                bubble: *gap_bubble,
                tip: *gap_tip,
                drop_repairable_bubbles: !*gap_keep_repairable,
            };
            opts.build_opts.k = *kmer_size;
            opts.build_opts.num_threads = *num_threads;
            opts.build_opts.prefix_precalc = *prefix_precalc;
            opts.build_opts.mem_gb = *mem_gb;
            opts.build_opts.temp_dir = Some(temp_dir.clone());

            info!("Choosing the index k-mers from {} inputs...", in_files.len());
            let read = |file: &String| read_fastx_file(file).into_iter().map(|(_, seq)| seq).collect::<Vec<Vec<u8>>>();
            let result = ntcomp::tigs::build(&in_files, &read, &opts, &|msg| info!("{}", msg));
            if own_temp_dir {
                let _ = std::fs::remove_dir_all(&temp_dir);
            }
            let bases: usize = result.tigs.iter().map(|t| t.len()).sum();
            info!("{} tigs, {} bases, {} canonical k-mers ({} in the core); {} tigs reversed by the layout",
                  result.tigs.len(), bases, result.stats.kmers, result.stats.core_kmers, result.stats.reversed);

            let tigs_path = format!("{}.tigs.tac", output_prefix);
            info!("Writing {}...", tigs_path);
            let out = BufWriter::new(File::create(&tigs_path).unwrap_or_else(|_| panic!("Could not create {}", tigs_path)));
            ntcomp::tigs::write_tigs(&result.tigs, *kmer_size, out).unwrap_or_else(|e| panic!("Could not write {}: {}", tigs_path, e)).flush().unwrap();

            if *fasta {
                let fasta_path = format!("{}.tigs.fasta", output_prefix);
                info!("Writing {}...", fasta_path);
                let mut out = BufWriter::new(File::create(&fasta_path).unwrap_or_else(|_| panic!("Could not create {}", fasta_path)));
                for (i, tig) in result.tigs.iter().enumerate() {
                    writeln!(out, ">{}", i).unwrap();
                    out.write_all(tig).unwrap();
                    out.write_all(b"\n").unwrap();
                }
            }
            if *index {
                info!("Serializing the index to {}.sbwt and {}.lcs ...", output_prefix, output_prefix);
                kbo::index::serialize_sbwt(output_prefix, &result.sbwt, &result.lcs);
            }
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

            let _ = ntcomp::nlz::legacy::write_header(&mut stdout);

            info!("Encoding fastX data...");
            let mut reader = needletail::parse_fastx_file(query_file).unwrap_or_else(|_| panic!("Expected valid fastX file"));
            let mut dictionaries: Vec<Vec<(usize, std::ops::Range<usize>)>> = Vec::new();
            let mut num_records = 0;

            while let Some(rec) = read_from_fastx_parser(&mut *reader) {
                let seqrec = rec.normalize(true);
                num_records += 1;

                dictionaries.push(ntcomp::nlz::parse(&seqrec, &sbwt, &lcs).unwrap().0);

                if num_records % block_size == 0 {
                    let u64_encoding = dictionaries.iter().flat_map(|x| ntcomp::nlz::legacy::encode::encode_dictionary(x, &sbwt).unwrap()).collect::<Vec<u64>>();
                    let _ = ntcomp::nlz::legacy::write_block_to(&u64_encoding, block_size, &mut stdout);
                    dictionaries.clear();
                }
            }
            if num_records % block_size != 0 {
                let u64_encoding = dictionaries.iter().flat_map(|x| ntcomp::nlz::legacy::encode::encode_dictionary(x, &sbwt).unwrap()).collect::<Vec<u64>>();
                let _ = ntcomp::nlz::legacy::write_block_to(&u64_encoding, num_records % block_size, &mut stdout);
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
            let _ = ntcomp::nlz::legacy::read_header(&mut conn);

            info!("Decoding encoded data...");
            let mut i = 0;
            while let Ok(encoding) = ntcomp::nlz::legacy::read_block(&mut conn) {
                let records = ntcomp::nlz::legacy::decode_sequence(&encoding, &sbwt);
                records.iter().for_each(|nucleotides| {
                    let _ = writeln!(&mut stdout, ">seq.{}", i + 1);
                    let _ = writeln!(&mut stdout,
                                     "{}", nucleotides.iter().map(|x| *x as char).collect::<String>());
                    let _ = stdout.flush();
                    i += 1;
                });
            }
        },
        Some(cli::Commands::Graph {
            query_files,
            index_prefix,
            decompress,
            fasta,
            fasta_columns,
            collapse_chains,
            reorient,
            reorient_sampling,
            edit_snps,
            edit_snps_max_shared,
            num_threads,
            batch_size,
            colors,
            contigs,
        }) => {
            init_log(2);
            let mut stdout = BufWriter::new(std::io::stdout());
            let batch_size = *batch_size as usize;

            rayon::ThreadPoolBuilder::new()
                .num_threads(*num_threads)
                .build_global()
                .unwrap();

            let (sbwt, lcs) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());

            if *decompress {
                assert!(query_files.len() == 1);
                let mut input = std::io::BufReader::new(File::open(&query_files[0]).unwrap());

                let decoder = ntcomp::graph::decoder::GraphDecoder::new(&mut input).unwrap();
                let header = decoder.header();
                let (graph, choices, hashes, contig_names, file_ranges, reversed, snp_edits) = decoder.decode_from_read_with_index(&sbwt, Some(&lcs)).unwrap();

                // Contigs to decode: all by default, otherwise the requested ones in the given order.
                let mut want_colors: Vec<u32> = colors.clone();
                if !contigs.is_empty() && contig_names.is_none() {
                    eprintln!("ntcomp: archive {} has no contig names", query_files[0].to_string_lossy());
                    std::process::exit(1);
                }
                for name in contigs {
                    let matches: Vec<u32> = contig_names.iter().flatten().enumerate()
                        .filter(|(_, x)| x.as_slice() == name.as_bytes() || x.split(|c| *c == b' ').next() == Some(name.as_bytes()))
                        .map(|(i, _)| i as u32)
                        .collect();
                    if matches.is_empty() {
                        eprintln!("ntcomp: contig {} not found in archive {}", name, query_files[0].to_string_lossy());
                        std::process::exit(1);
                    }
                    want_colors.extend(matches);
                }
                if want_colors.is_empty() {
                    want_colors = (0..header.n_queries).collect();
                }
                if let Some(bad) = want_colors.iter().find(|&&c| c >= header.n_queries) {
                    eprintln!("ntcomp: color {} out of range (archive {} contains {} contigs)", bad, query_files[0].to_string_lossy(), header.n_queries);
                    std::process::exit(1);
                }

                if !fasta {
                    eprintln!("color\taccession\tcontig\tdecoded_len\thash_match\tseconds");
                }
                // Colours are decoded in parallel in chunks of about `batch_size` input
                // files' worth of contigs, and written in the requested order.
                let n_files = file_ranges.as_ref().map_or(0, |ranges| ranges.len()).max(1);
                let chunk_size = batch_size * (header.n_queries as usize).div_ceil(n_files).max(1);
                // Repeated SNPs are located by contexts from the edited sequences of
                // earlier colours: the positions needed, by colour.
                let mut wanted_contexts: std::collections::HashMap<u32, Vec<u32>> = std::collections::HashMap::new();
                if let Some(edits) = &snp_edits {
                    for &color in &want_colors {
                        for (f, pos) in edits.context_sources(color) { wanted_contexts.entry(f).or_default().push(pos); }
                    }
                }
                let mut contexts: SnpContexts = SnpContexts::new();
                let record = |contexts: &mut SnpContexts, color: u32, edited: &[u8]| {
                    for &pos in wanted_contexts.get(&color).map_or(&[][..], |v| v.as_slice()) {
                        if let Some(context) = ntcomp::graph::snps::context_at(edited, pos as usize) { contexts.insert((color, pos), context); }
                    }
                };
                for chunk in want_colors.chunks(chunk_size) {
                    let start_time = std::time::Instant::now();
                    let edited: Vec<Result<Vec<u8>, String>> = chunk.par_iter().map(|&color| {
                        decode_edited(color, &graph, &choices, &sbwt, &lcs)
                    }).collect();
                    for (&color, sequence) in chunk.iter().zip(&edited) {
                        if let Ok(sequence) = sequence { record(&mut contexts, color, sequence); }
                    }
                    // Contexts from colours not decoded yet (random access).
                    if let Some(edits) = &snp_edits {
                        let mut missing: Vec<u32> = chunk.iter().flat_map(|&c| edits.context_sources(c)).filter(|key| !contexts.contains_key(key)).map(|(f, _)| f).collect();
                        missing.sort_unstable();
                        missing.dedup();
                        let sources: Vec<(u32, Result<Vec<u8>, String>)> = missing.par_iter().map(|&f| (f, decode_edited(f, &graph, &choices, &sbwt, &lcs))).collect();
                        for (f, sequence) in sources {
                            let sequence = sequence.unwrap_or_else(|e| panic!("decoding color {f} failed: {e}"));
                            record(&mut contexts, f, &sequence);
                        }
                    }
                    let seconds = start_time.elapsed().as_secs_f64() / chunk.len() as f64;
                    let decoded: Vec<Result<(Vec<u8>, Option<bool>, f64), String>> = chunk.par_iter().zip(edited.into_par_iter()).map(|(&color, sequence)| {
                        finish_color(color, sequence?, &hashes, &reversed, &snp_edits, &contexts).map(|(s, h)| (s, h, seconds))
                    }).collect();

                    for (&color, decoded) in chunk.iter().zip(decoded) {
                        let idx = color as usize;
                        let (sequence, hash_match, seconds) = decoded.unwrap_or_else(|e| panic!("decoding color {color} failed: {e}"));
                        assert!(hash_match != Some(false), "hash mismatch for color {color}");

                        let contig_name: Vec<u8> = contig_names.as_ref().map_or_else(|| format!("color_{}", color).into_bytes(), |names| names[idx].clone());
                        if !fasta {
                            let file_name = file_ranges.as_ref()
                                .and_then(|ranges| ranges.iter().position(|(_, range)| range.contains(&color)).map(|i| String::from_utf8_lossy(&ranges[i].0).to_string()))
                                .unwrap_or("NA".to_string());
                            let hash_match = hash_match.map_or("NA".to_string(), |x| x.to_string());
                            eprintln!("{}\t{}\t{}\t{}\t{}\t{:.3}", color, file_name, String::from_utf8_lossy(&contig_name), sequence.len(), hash_match, seconds);
                        } else {
                            stdout.write_all(b">").unwrap();
                            stdout.write_all(&contig_name).unwrap();
                            stdout.write_all(b"\n").unwrap();
                            sequence.chunks(*fasta_columns).for_each(|chunk| {
                                stdout.write_all(chunk).unwrap();
                                stdout.write_all(b"\n").unwrap();
                            });
                        }
                    }
                }
                stdout.flush().unwrap();
            } else {
                let n_queries = query_files.len();

                let progress = ProgressBar::new(n_queries as u64);
                progress.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}").unwrap());

                let mut builder = ntcomp::graph::GraphBuilder::new()
                    .collapse_chains(*collapse_chains)
                    .reorient_contigs(*reorient)
                    .reorient_sampling(*reorient_sampling)
                    .edit_snps(*edit_snps)
                    .edit_snps_max_shared(edit_snps_max_shared.map_or(ntcomp::graph::snps::DEFAULT_MAX_SHARED, |max| max as usize));

                // First pass for SNP editing: choose the strands and call the SNPs of every
                // contig, `batch_size` input files at a time.
                if *edit_snps {
                    let progress = ProgressBar::new(n_queries as u64);
                    progress.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} calling SNPs").unwrap());
                    for batch in query_files.chunks(batch_size) {
                        let contigs: Vec<Vec<(Vec<u8>, Vec<u8>)>> = batch.par_iter().map(read_contigs).collect();
                        let sequences: Vec<&[u8]> = contigs.iter().flatten().map(|(_, seq)| seq.as_slice()).collect();
                        let called = builder.call_edits_batch(&sequences, &sbwt, &lcs);
                        info!("Called SNPs of colors {}..={}", called.first().map_or(0, |c| *c), called.last().map_or(0, |c| *c));
                        progress.inc(batch.len() as u64);
                    }
                    progress.finish();
                }

                // Hash of every (normalized) input contig, checked after decoding
                let mut hashes: Vec<ContigHash> = Vec::new();

                let mut contig_names: Vec<Vec<u8>> = Vec::with_capacity(n_queries);
                let mut file_ranges: Vec<(Vec<u8>, core::ops::Range<u32>)> = Vec::with_capacity(n_queries);

                let mut color: u32 = 0;
                for batch in query_files.chunks(batch_size) {
                    let contigs: Vec<Vec<(Vec<u8>, Vec<u8>)>> = batch.par_iter().map(read_contigs).collect();
                    let sequences: Vec<&[u8]> = contigs.iter().flatten().map(|(_, seq)| seq.as_slice()).collect();
                    hashes.par_extend(sequences.par_iter().map(|seq| contig_hash(seq)));

                    let added = builder.add_contigs(&sequences, &sbwt, &lcs).unwrap();
                    for &added in &added {
                        info!("Encoded color {}{}", added, if builder.reversed()[added as usize] { " (reverse complemented)" } else { "" });
                    }

                    for (query_file, file_contigs) in batch.iter().zip(&contigs) {
                        let start = color;
                        for (name, _) in file_contigs {
                            contig_names.push(name.clone());
                            assert_eq!(added[(color - added[0]) as usize], color);
                            color += 1;
                        }
                        let filename: Vec<u8> = query_file.file_prefix().unwrap().as_encoded_bytes().to_vec();
                        file_ranges.push((filename, start..color));
                    }

                    progress.inc(batch.len() as u64);
                }
                progress.finish();

                let reversed = builder.reversed().to_vec();
                let snp_edits = builder.snp_edits();
                let (graph, choices) = builder.finish();

                let mut encoder = ntcomp::graph::encoder::GraphEncoder::new(graph, choices)
                    .with_hashes(hashes)
                    .with_contig_names(contig_names)
                    .with_file_ranges(file_ranges)
                    .with_reversed(reversed)
                    .with_index(&sbwt, &lcs);
                if let Some(snp_edits) = snp_edits {
                    encoder = encoder.with_snp_edits(snp_edits);
                }

                let sizes = encoder.encode_to_write(&mut stdout).unwrap();
                info!("{}", sizes);
            }

            let _ = stdout.flush();
        },
        #[cfg(feature = "graph-view")]
        Some(cli::Commands::GraphView {
            input_path,
            index_prefix,
            color,
        }) => {
            init_log(2);

            let (sbwt, lcs) = kbo::index::load_sbwt(index_prefix.as_ref().unwrap());

            let mut input = std::io::BufReader::new(File::open(&input_path).unwrap());
            let decoder = ntcomp::graph::decoder::GraphDecoder::new(&mut input).unwrap();
            let header = decoder.header();
            let (stored_graph, choices, _hashes, _contig_names, _file_ranges, _reversed, _snp_edits) = decoder.decode_from_read_with_index(&sbwt, Some(&lcs)).unwrap();

            // Colour sets are not stored, recover them by walking every colour.
            let graph: petgraph::graph::Graph<u32, ntcomp::graph::ColexGraphEdge> = {
                let csr = stored_graph.into_csr();
                let colorsets = csr.colorsets(&choices, header.n_queries).unwrap();
                csr.to_petgraph(&colorsets)
            };

            let colors: Vec<u32> = (0..header.n_queries).collect();

            let native_options = eframe::NativeOptions::default();
            let ret = eframe::run_native(
                "Interactive Graph Viewer",
                native_options,
                Box::new(|_| Ok(Box::new
                                 (ntcomp::graphview::BasicApp::new(
                                     graph,
                                     colors[*color],
                                 )))),
            );
            ret.unwrap();
        },
        None => {},
    }
}
