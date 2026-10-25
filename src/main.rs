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

use cli::{BuildArgs, Commands, DecodeArgs, EncodeArgs, GraphArgs, TigsArgs};
#[cfg(feature = "graph-view")]
use cli::GraphViewArgs;

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

/// `ntcomp build`: an index of FASTA, FASTQ or tig files.
fn build(
    args: &BuildArgs,
) {
    init_log(if args.verbose { 2 } else { 1 });

    let in_files = args.input_files().or_exit();
    info!("Building SBWT index from {} files...", in_files.len());
    let mut seq_data: Vec<Vec<u8>> = Vec::new();
    for file in &in_files {
        seq_data.append(&mut ntcomp::tigs::read_sequences(file).or_exit());
    }

    let index = Index::build(&seq_data, &args.build_opts());

    info!("Serializing SBWT index to {}.sbwt ...", args.output_prefix);
    info!("Serializing LCS array to {}.lcs ...", args.output_prefix);
    index.save(&args.output_prefix).or_exit();
}

/// `ntcomp tigs`: chooses the index k-mers of a collection and writes them
/// as tigs.
fn tigs(
    args: &TigsArgs,
) {
    init_log(if args.verbose { 2 } else { 1 });

    let in_files = args.input_files().or_exit();
    let (temp_dir, own_temp_dir) = args.temp_dir();
    std::fs::create_dir_all(&temp_dir).map_err(|e| format!("creating {}: {}", temp_dir, e)).or_exit();
    let opts = args.tig_opts(&temp_dir);

    info!("Choosing the index k-mers from {} inputs...", in_files.len());
    let read = |file: &String| fastx::read_contigs(file).or_exit().into_iter().map(|(_, sequence)| sequence).collect::<Vec<Vec<u8>>>();
    let result = ntcomp::tigs::build(&in_files, &read, &opts, &|msg| info!("{}", msg));
    if own_temp_dir {
        let _ = std::fs::remove_dir_all(&temp_dir);
    }
    let bases: usize = result.tigs.iter().map(|t| t.len()).sum();
    info!("{} tigs, {} bases, {} canonical k-mers ({} in the core); {} tigs reversed by the layout",
          result.tigs.len(), bases, result.stats.kmers, result.stats.core_kmers, result.stats.reversed);

    let output_prefix = &args.output_prefix;
    let tigs_path = format!("{}.tigs.tac", output_prefix);
    info!("Writing {}...", tigs_path);
    let out = BufWriter::new(File::create(&tigs_path).map_err(|e| format!("creating {}: {}", tigs_path, e)).or_exit());
    ntcomp::tigs::write_tigs(&result.tigs, args.kmer_size, out).map_err(|e| format!("writing {}: {}", tigs_path, e)).or_exit().flush().or_exit();

    if args.fasta {
        let fasta_path = format!("{}.tigs.fasta", output_prefix);
        info!("Writing {}...", fasta_path);
        let mut out = BufWriter::new(File::create(&fasta_path).map_err(|e| format!("creating {}: {}", fasta_path, e)).or_exit());
        for (i, tig) in result.tigs.iter().enumerate() {
            fastx::write_fasta(&mut out, i.to_string().as_bytes(), tig, usize::MAX).or_exit();
        }
        out.flush().or_exit();
    }
    if args.index {
        info!("Serializing the index to {}.sbwt and {}.lcs ...", output_prefix, output_prefix);
        result.index.save(output_prefix).or_exit();
    }
}

/// `ntcomp encode`: a FASTA or FASTQ file in the plain NLZ format.
fn encode(
    args: &EncodeArgs,
) {
    init_log(2);
    let mut stdout = BufWriter::new(std::io::stdout());
    info!("Loading SBWT index...");
    let index = Index::load(&args.index_prefix).or_exit();

    info!("Encoding fastX data...");
    let contigs = fastx::read_contigs(&args.query_file).or_exit();
    ntcomp::nlz::legacy::write_file(contigs.iter().map(|(_, sequence)| sequence), &index, &mut stdout).or_exit();
    stdout.flush().or_exit();
}

/// `ntcomp decode`: decodes a plain NLZ file.
fn decode(
    args: &DecodeArgs,
) {
    init_log(2);
    let mut stdout = BufWriter::new(std::io::stdout());
    let index = Index::load(&args.index_prefix).or_exit();
    let input_path = &args.input_path;
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
}

/// `ntcomp graph`: compresses assemblies into an archive, or decompresses
/// one (`-d`).
fn graph(
    args: &GraphArgs,
) {
    init_log(2);
    let mut stdout = BufWriter::new(std::io::stdout());

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.num_threads)
        .build_global()
        .unwrap();

    let index = Index::load(&args.index_prefix).or_exit();
    if args.decompress {
        graph_decode(args, &index, &mut stdout);
    } else {
        graph_encode(args, &index, &mut stdout);
    }
    stdout.flush().or_exit();
}

/// Writes the archive of `args.query_files` to `out`.
fn graph_encode<W: Write>(
    args: &GraphArgs,
    index: &Index,
    out: &mut W,
) {
    let (query_files, batch_size) = (&args.query_files, args.batch_size as usize);
    let n_queries = query_files.len();

    let progress = ProgressBar::new(n_queries as u64);
    progress.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}").unwrap());

    let mut writer = ArchiveWriter::new(args.archive_opts(), index);

    for batch in query_files.chunks(batch_size) {
        let added = writer.add_files(&read_input_files(batch)).or_exit();
        for &added in &added {
            info!("Encoded color {}{}", added, if writer.reversed()[added as usize] { " (reverse complemented)" } else { "" });
        }
        progress.inc(batch.len() as u64);
    }
    progress.finish();

    let sizes = writer.write(out).or_exit();
    info!("{}", sizes);
}

/// Decodes the requested contigs of the archive `args.query_files[0]`
/// (all by default) to `out` as FASTA (`--fasta`), or reports them; exits
/// with 1 after all of them if one failed or did not match its hash.
fn graph_decode<W: Write>(
    args: &GraphArgs,
    index: &Index,
    out: &mut W,
) {
    if args.query_files.len() != 1 {
        eprintln!("ntcomp: decompressing takes one archive, got {}", args.query_files.len());
        std::process::exit(1);
    }
    let archive_path = args.query_files[0].to_string_lossy();
    let in_archive = |e: &dyn std::fmt::Display| format!("archive {}: {}", archive_path, e);
    let mut input = std::io::BufReader::new(File::open(&args.query_files[0]).map_err(|e| in_archive(&e)).or_exit());
    let archive = ArchiveReader::open(&mut input, index).map_err(|e| in_archive(&e)).or_exit();

    // Contigs to decode: all by default, otherwise the requested ones in the given order.
    let mut want_colors: Vec<u32> = args.colors.clone();
    for name in &args.contigs {
        want_colors.extend(archive.find_contig(name).map_err(|e| in_archive(&e)).or_exit());
    }
    if want_colors.is_empty() {
        want_colors = (0..archive.n_contigs()).collect();
    }
    let mut decoder = archive.decoder(&want_colors).map_err(|e| in_archive(&e)).or_exit();

    if !args.fasta {
        eprintln!("color\taccession\tcontig\tdecoded_len\thash_match\tseconds");
    }
    // Colours are decoded in parallel in chunks, and written in the requested order.
    let mut failed = false;
    for chunk in want_colors.chunks(args.chunk_size(archive.n_contigs(), archive.n_files())) {
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
            if !args.fasta {
                let file_name = archive.file_name(color).map_or("NA".to_string(), |name| String::from_utf8_lossy(name).to_string());
                let hash_match = contig.hash_match.map_or("NA".to_string(), |x| x.to_string());
                eprintln!("{}\t{}\t{}\t{}\t{}\t{:.3}", color, file_name, String::from_utf8_lossy(&contig_name), contig.sequence.len(), hash_match, seconds);
            } else {
                fastx::write_fasta(out, &contig_name, &contig.sequence, args.fasta_columns).or_exit();
            }
        }
    }
    if failed {
        out.flush().or_exit();
        std::process::exit(1);
    }
}

/// `ntcomp graph-view`: shows the graph of an archive.
#[cfg(feature = "graph-view")]
fn graph_view(
    args: &GraphViewArgs,
) {
    init_log(2);

    let index = Index::load(&args.index_prefix).or_exit();
    let input_path = &args.input_path;
    let mut input = std::io::BufReader::new(File::open(input_path).map_err(|e| format!("{}: {}", input_path.display(), e)).or_exit());
    let archive = ArchiveReader::open(&mut input, &index).or_exit();
    if args.color >= archive.n_contigs() as usize {
        eprintln!("ntcomp: color {} out of range ({} contigs)", args.color, archive.n_contigs());
        std::process::exit(1);
    }
    // Colour sets are not stored, recover them by walking every colour.
    let (csr, colorsets) = archive.into_colored_graph().or_exit();
    let graph = ntcomp::graphview::to_petgraph(csr, &colorsets);

    let native_options = eframe::NativeOptions::default();
    let ret = eframe::run_native(
        "Interactive Graph Viewer",
        native_options,
        Box::new(|_| Ok(Box::new
                         (ntcomp::graphview::BasicApp::new(
                             graph,
                             args.color as u32,
                         )))),
    );
    ret.unwrap();
}

fn main() {
    match &cli::Cli::parse().command {
        Some(Commands::Build(args)) => build(args),
        Some(Commands::Tigs(args)) => tigs(args),
        Some(Commands::Encode(args)) => encode(args),
        Some(Commands::Decode(args)) => decode(args),
        Some(Commands::Graph(args)) => graph(args),
        #[cfg(feature = "graph-view")]
        Some(Commands::GraphView(args)) => graph_view(args),
        None => {},
    }
}
