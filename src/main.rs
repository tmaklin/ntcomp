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
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

use clap::Parser;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use log::info;
use rayon::prelude::*;

use ntcomp::archive::InputFile;
use ntcomp::fastx;
use ntcomp::{ArchiveReader, ArchiveWriter, Index};

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

/// Exits with an error message instead of a panic.
trait OrExit<T> {
    fn or_exit(self) -> T;
}

impl<T, D: std::fmt::Display> OrExit<T> for Result<T, D> {
    fn or_exit(self) -> T {
        self.unwrap_or_else(|e| {
            eprintln!("ntcomp: {}", e);
            std::process::exit(1)
        })
    }
}

/// The input files: `files`, then those listed in `input_list`.
fn input_files(
    files: &[String],
    input_list: &Option<String>,
) -> Vec<String> {
    let mut in_files = files.to_vec();
    if let Some(list) = input_list {
        in_files.extend(cli::read_input_list(list, b'\t').or_exit());
    }
    in_files
}

/// Reads `files` in parallel, as archive inputs named by their file name
/// without extensions.
fn read_input_files(
    files: &[std::path::PathBuf],
) -> Vec<InputFile> {
    files.par_iter().map(|path| {
        let contigs = fastx::read_contigs(path).map_err(|e| e.to_string())?;
        let name = path.file_prefix().unwrap_or_default().as_encoded_bytes().to_vec();
        Ok((name, contigs))
    }).collect::<Result<Vec<InputFile>, String>>().or_exit()
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

            let build_opts = ntcomp::index::BuildOpts {
                k: *kmer_size,
                num_threads: *num_threads,
                prefix_precalc: *prefix_precalc,
                mem_gb: *mem_gb,
                dedup_batches: *dedup_batches,
                temp_dir: temp_dir.clone(),
            };

            let in_files = input_files(seq_files, input_list);
            info!("Building SBWT index from {} files...", in_files.len());
            let mut seq_data: Vec<Vec<u8>> = Vec::new();
            for file in &in_files {
                seq_data.append(&mut fastx::read_sequences(file).or_exit());
            }

            let index = Index::build(&seq_data, &build_opts);

            let output_prefix = output_prefix.as_ref().unwrap();
            info!("Serializing SBWT index to {}.sbwt ...", output_prefix);
            info!("Serializing LCS array to {}.lcs ...", output_prefix);
            index.save(output_prefix).or_exit();
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

            let in_files = input_files(seq_files, input_list);

            // The indexes are built on disk: in memory, sorting the k-mers takes
            // k / 4 bytes each.
            let (temp_dir, own_temp_dir) = match temp_dir {
                Some(dir) => (dir.clone(), false),
                None => (format!("{}.tmp", output_prefix), true),
            };
            std::fs::create_dir_all(&temp_dir).map_err(|e| format!("creating {}: {}", temp_dir, e)).or_exit();

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
            let read = |file: &String| fastx::read_contigs(file).or_exit().into_iter().map(|(_, sequence)| sequence).collect::<Vec<Vec<u8>>>();
            let result = ntcomp::tigs::build(&in_files, &read, &opts, &|msg| info!("{}", msg));
            if own_temp_dir {
                let _ = std::fs::remove_dir_all(&temp_dir);
            }
            let bases: usize = result.tigs.iter().map(|t| t.len()).sum();
            info!("{} tigs, {} bases, {} canonical k-mers ({} in the core); {} tigs reversed by the layout",
                  result.tigs.len(), bases, result.stats.kmers, result.stats.core_kmers, result.stats.reversed);

            let tigs_path = format!("{}.tigs.tac", output_prefix);
            info!("Writing {}...", tigs_path);
            let out = BufWriter::new(File::create(&tigs_path).map_err(|e| format!("creating {}: {}", tigs_path, e)).or_exit());
            ntcomp::tigs::write_tigs(&result.tigs, *kmer_size, out).map_err(|e| format!("writing {}: {}", tigs_path, e)).or_exit().flush().or_exit();

            if *fasta {
                let fasta_path = format!("{}.tigs.fasta", output_prefix);
                info!("Writing {}...", fasta_path);
                let mut out = BufWriter::new(File::create(&fasta_path).map_err(|e| format!("creating {}: {}", fasta_path, e)).or_exit());
                for (i, tig) in result.tigs.iter().enumerate() {
                    fastx::write_fasta(&mut out, i.to_string().as_bytes(), tig, usize::MAX).or_exit();
                }
                out.flush().or_exit();
            }
            if *index {
                info!("Serializing the index to {}.sbwt and {}.lcs ...", output_prefix, output_prefix);
                result.index.save(output_prefix).or_exit();
            }
        },
        Some(cli::Commands::Encode {
            query_file,
            index_prefix,
        }) => {
            init_log(2);
            let mut stdout = BufWriter::new(std::io::stdout());
            info!("Loading SBWT index...");
            let index = Index::load(index_prefix.as_ref().unwrap()).or_exit();

            info!("Encoding fastX data...");
            let contigs = fastx::read_contigs(query_file).or_exit();
            ntcomp::nlz::legacy::write_file(contigs.iter().map(|(_, sequence)| sequence), &index, &mut stdout).or_exit();
            stdout.flush().or_exit();
        },

        Some(cli::Commands::Decode {
            input_path,
            index_prefix,
        }) => {
            init_log(2);
            let mut stdout = BufWriter::new(std::io::stdout());
            let index = Index::load(index_prefix.as_ref().unwrap()).or_exit();
            let mut input = File::open(input_path).map_err(|e| format!("{}: {}", input_path.display(), e)).or_exit();

            info!("Decoding encoded data...");
            let mut i = 0;
            ntcomp::nlz::legacy::read_file(&mut input, &index, |nucleotides| {
                i += 1;
                writeln!(&mut stdout, ">seq.{}", i)?;
                stdout.write_all(&nucleotides)?;
                stdout.write_all(b"\n")?;
                Ok(())
            }).or_exit();
            stdout.flush().or_exit();
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

            let index = Index::load(index_prefix.as_ref().unwrap()).or_exit();

            if *decompress {
                if query_files.len() != 1 {
                    eprintln!("ntcomp: decompressing takes one archive, got {}", query_files.len());
                    std::process::exit(1);
                }
                let archive_path = query_files[0].to_string_lossy();
                let in_archive = |e: &dyn std::fmt::Display| format!("archive {}: {}", archive_path, e);
                let mut input = std::io::BufReader::new(File::open(&query_files[0]).map_err(|e| in_archive(&e)).or_exit());
                let archive = ArchiveReader::open(&mut input, &index).map_err(|e| in_archive(&e)).or_exit();

                // Contigs to decode: all by default, otherwise the requested ones in the given order.
                let mut want_colors: Vec<u32> = colors.clone();
                for name in contigs {
                    want_colors.extend(archive.find_contig(name).map_err(|e| in_archive(&e)).or_exit());
                }
                if want_colors.is_empty() {
                    want_colors = (0..archive.n_contigs()).collect();
                }
                let mut decoder = archive.decoder(&want_colors).map_err(|e| in_archive(&e)).or_exit();

                if !fasta {
                    eprintln!("color\taccession\tcontig\tdecoded_len\thash_match\tseconds");
                }
                // Colours are decoded in parallel in chunks of about `batch_size` input
                // files' worth of contigs, and written in the requested order.
                let chunk_size = batch_size * (archive.n_contigs() as usize).div_ceil(archive.n_files().max(1)).max(1);
                let mut failed = false;
                for chunk in want_colors.chunks(chunk_size) {
                    let start_time = std::time::Instant::now();
                    let decoded = decoder.decode(chunk).map_err(|e| in_archive(&e)).or_exit();
                    let seconds = start_time.elapsed().as_secs_f64() / chunk.len() as f64;

                    for (&color, decoded) in chunk.iter().zip(decoded) {
                        let contig = match decoded {
                            Ok(contig) => contig,
                            Err(e) => {
                                eprintln!("ntcomp: {}", in_archive(&e));
                                failed = true;
                                continue
                            },
                        };
                        failed |= contig.hash_match == Some(false);
                        let contig_name: Vec<u8> = archive.contig_name(color).map_or_else(|| format!("color_{}", color).into_bytes(), |name| name.to_vec());
                        if !fasta {
                            let file_name = archive.file_name(color).map_or("NA".to_string(), |name| String::from_utf8_lossy(name).to_string());
                            let hash_match = contig.hash_match.map_or("NA".to_string(), |x| x.to_string());
                            eprintln!("{}\t{}\t{}\t{}\t{}\t{:.3}", color, file_name, String::from_utf8_lossy(&contig_name), contig.sequence.len(), hash_match, seconds);
                        } else {
                            fastx::write_fasta(&mut stdout, &contig_name, &contig.sequence, *fasta_columns).or_exit();
                        }
                    }
                }
                stdout.flush().or_exit();
                if failed {
                    std::process::exit(1);
                }
            } else {
                let n_queries = query_files.len();

                let progress = ProgressBar::new(n_queries as u64);
                progress.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}").unwrap());

                let builder = ntcomp::graph::GraphBuilder::new()
                    .collapse_chains(*collapse_chains)
                    .reorient_contigs(*reorient)
                    .reorient_sampling(*reorient_sampling)
                    .edit_snps(*edit_snps)
                    .edit_snps_max_shared(edit_snps_max_shared.map_or(ntcomp::graph::snps::DEFAULT_MAX_SHARED, |max| max as usize));
                let mut writer = ArchiveWriter::new(builder, &index);

                // First pass for SNP editing: choose the strands and call the SNPs of every
                // contig, `batch_size` input files at a time.
                if *edit_snps {
                    let progress = ProgressBar::new(n_queries as u64);
                    progress.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} calling SNPs").unwrap());
                    for batch in query_files.chunks(batch_size) {
                        let called = writer.call_edits(&read_input_files(batch));
                        info!("Called SNPs of colors {}..={}", called.first().map_or(0, |c| *c), called.last().map_or(0, |c| *c));
                        progress.inc(batch.len() as u64);
                    }
                    progress.finish();
                }

                for batch in query_files.chunks(batch_size) {
                    let added = writer.add_files(&read_input_files(batch)).or_exit();
                    for &added in &added {
                        info!("Encoded color {}{}", added, if writer.reversed()[added as usize] { " (reverse complemented)" } else { "" });
                    }
                    progress.inc(batch.len() as u64);
                }
                progress.finish();

                let sizes = writer.write(&mut stdout).or_exit();
                info!("{}", sizes);
            }

            stdout.flush().or_exit();
        },
        #[cfg(feature = "graph-view")]
        Some(cli::Commands::GraphView {
            input_path,
            index_prefix,
            color,
        }) => {
            init_log(2);

            let index = Index::load(index_prefix.as_ref().unwrap()).or_exit();

            let mut input = std::io::BufReader::new(File::open(input_path).map_err(|e| format!("{}: {}", input_path.display(), e)).or_exit());
            let archive = ArchiveReader::open(&mut input, &index).or_exit();
            if *color >= archive.n_contigs() as usize {
                eprintln!("ntcomp: color {} out of range ({} contigs)", color, archive.n_contigs());
                std::process::exit(1);
            }
            // Colour sets are not stored, recover them by walking every colour.
            let graph = archive.to_petgraph().or_exit();

            let native_options = eframe::NativeOptions::default();
            let ret = eframe::run_native(
                "Interactive Graph Viewer",
                native_options,
                Box::new(|_| Ok(Box::new
                                 (ntcomp::graphview::BasicApp::new(
                                     graph,
                                     *color as u32,
                                 )))),
            );
            ret.unwrap();
        },
        None => {},
    }
}
