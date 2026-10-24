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
use std::path::PathBuf;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(version)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    // Build the compression dictionary
    Build {
        // Input fasta or fastq sequence file(s)
        #[arg(group = "input", required = true, help = "Sequence data file(s).")]
        seq_files: Vec<String>,

        #[arg(short = 'l', long = "input-list", group = "input", required = true, help_heading = "Input", help = "File with paths or tab separated name and path on each line.")]
        input_list: Option<String>,

        // Outputs
        #[arg(short = 'o', long = "output-prefix", required = true, help_heading = "Output", help = "Prefix for output files <prefix>.sbwt and <prefix>.lcs.")]
        output_prefix: Option<String>,

        // Build parameters
        // // k-mer size
        #[arg(short = 'k', default_value_t = 31, help_heading = "Build options", help = "k-mer size, larger values are slower and use more space.")]
        kmer_size: usize,
        // // prefix precalc
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options", help = "Length of precalculated prefixes included in the index.")]
        prefix_precalc: usize,
        // // deduplicate k-mer batches
        #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options", help = "Deduplicate k-mer batches to save some memory.")]
        dedup_batches: bool,

        // Resources
        // // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1)]
        num_threads: usize,
        // // Memory in GB
        #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options", help = "Memory available when building on temp disk space (in gigabytes).")]
        mem_gb: usize,
        // // Temporary directory
        #[arg(long = "temp-dir", required = false, help_heading = "Build options", help = "Build on temporary disk space at this path instead of in-memory.")]
        temp_dir: Option<String>,

        // Verbosity
        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },

    // Choose the index k-mers from a collection and write them as tigs
    Tigs {
        // Input fasta or fastq sequence file(s)
        #[arg(group = "input", required = true, help = "Assembly file(s), one per input; the tigs are seeded in this order.")]
        seq_files: Vec<String>,

        #[arg(short = 'l', long = "input-list", group = "input", required = true, help_heading = "Input", help = "File with paths or tab separated name and path on each line.")]
        input_list: Option<String>,

        // Outputs
        #[arg(short = 'o', long = "output-prefix", required = true, help_heading = "Output", help = "Prefix for <prefix>.tigs.tac (and the outputs below).")]
        output_prefix: String,

        #[arg(long = "fasta", default_value_t = false, help_heading = "Output", help = "Also write the tigs to <prefix>.tigs.fasta.")]
        fasta: bool,

        #[arg(long = "index", default_value_t = false, help_heading = "Output", help = "Also write the index built from the tigs to <prefix>.sbwt and <prefix>.lcs.")]
        index: bool,

        // Tig parameters
        #[arg(short = 'k', default_value_t = 255, help_heading = "Tig options", help = "k-mer size.")]
        kmer_size: usize,

        #[arg(short = 'T', long = "min-assemblies", default_value_t = 5, help_heading = "Tig options", help = "Core: canonical k-mers present in at least this many inputs.")]
        min_assemblies: u32,

        #[arg(long = "gap-join", default_value = "1-inf", value_parser = parse_gap_range, help_heading = "Tig options", help = "Lengths (in k-mers) of gaps between two core tig ends to fill: <min>-<max|inf> or none.")]
        gap_join: ntcomp::tigs::GapRange,

        #[arg(long = "gap-one-sided", default_value = "1-inf", value_parser = parse_gap_range, help_heading = "Tig options", help = "Lengths of gaps with exactly one flank a core tig end to fill.")]
        gap_one_sided: ntcomp::tigs::GapRange,

        #[arg(long = "gap-bubble", default_value = "350-inf", value_parser = parse_gap_range, help_heading = "Tig options", help = "Lengths of gaps with neither flank a core tig end to fill.")]
        gap_bubble: ntcomp::tigs::GapRange,

        #[arg(long = "gap-tip", default_value = "250-inf", value_parser = parse_gap_range, help_heading = "Tig options", help = "Lengths of runs at contig ends to fill.")]
        gap_tip: ntcomp::tigs::GapRange,

        #[arg(long = "gap-keep-repairable", default_value_t = false, help_heading = "Tig options", help = "Also fill the bubbles that SNP and indel editing repairs against the core (left out by default).")]
        gap_keep_repairable: bool,

        #[arg(long = "orient-sampling", default_value_t = 8, help_heading = "Tig options", help = "Sample one in this many 31-mers for orienting the tigs.")]
        orient_sampling: u64,

        #[arg(long = "partitions", default_value_t = 16, help_heading = "Tig options", help = "Count the core in this many hash partitions (more use less memory but read the inputs more times).")]
        partitions: u64,

        // Index build parameters
        #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options", help = "Length of precalculated prefixes included in the index.")]
        prefix_precalc: usize,

        #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options", help = "Memory for building the indexes on temporary disk space (in gigabytes).")]
        mem_gb: usize,

        #[arg(long = "temp-dir", required = false, help_heading = "Build options", help = "Temporary disk space for building the indexes [default: <prefix>.tmp].")]
        temp_dir: Option<String>,

        // Resources
        #[arg(short = 't', long = "threads", default_value_t = 1, help = "Number of threads.")]
        num_threads: usize,

        #[arg(long = "verbose", default_value_t = false)]
        verbose: bool,
    },

    // Encode fastX data using an SBWT index
    Encode {
        // Inputs
        // // Input fasta or fastq query file(s)
        #[arg(group = "input", required = true, help = "Query file with sequence data.")]
        query_file: PathBuf,

        // // Prebuilt index
        #[arg(short = 'i', long = "index", help_heading = "Input", required = true, help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs")]
        index_prefix: Option<String>,

    },

    // Encode a collection of fastX files using an SBWT index
    Graph {
        // Inputs
        // // Input fasta or fastq query file(s)
        #[arg(group = "input", required = true, help = "Querys file with sequence data.")]
        query_files: Vec<PathBuf>,

        // // Prebuilt index
        #[arg(short = 'i', long = "index", help_heading = "Input", required = true, help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs")]
        index_prefix: Option<String>,

        #[arg(long = "fasta", default_value_t = false, help = ".fasta output")]
        fasta: bool,

        #[arg(long = "fasta-cols", default_value_t = 80, help = ".fasta columns")]
        fasta_columns: usize,

        // // Collapse linear chains
        #[arg(long = "collapse-chains", default_value_t = false, help = "Store linear chains of nodes as single edges.")]
        collapse_chains: bool,

        // // Choose contig strands
        #[arg(long = "reorient", default_value_t = false, help = "Store every contig on the strand that shares the most with the contigs before it.")]
        reorient: bool,

        #[arg(long = "reorient-sampling", default_value_t = 32, value_parser = clap::value_parser!(u64).range(1..), help = "Sample one in this many 31-mers for choosing the strands.")]
        reorient_sampling: u64,

        // // Edit SNPs towards the index
        #[arg(long = "edit-snps", default_value_t = false, help = "Edit SNPs and small indels against the index before encoding (reads the input twice).")]
        edit_snps: bool,


        // Resources
        // // Threads
        #[arg(short = 't', long = "threads", default_value_t = 1, help = "Number of threads.")]
        num_threads: usize,

        // // Batch size
        #[arg(long = "batch-size", default_value_t = 64, value_parser = clap::value_parser!(u64).range(1..), help = "Read and process this many input files at a time (when decoding: this many files' worth of contigs).")]
        batch_size: u64,

        // Decompression
        // // Toggle decompression
        #[arg(short = 'd', long = "decompress", default_value_t = false, help = "Force decompression", group = "mode")]
        decompress: bool,
        // // Retrieve specific contig(s) when decompressing
        #[arg(long = "color", help_heading = "Decompression", help = "Decode by contig index(es) (0-based).")]
        colors: Vec<u32>,

        #[arg(long = "contig", help_heading = "Decompression", help = "Decode by contig name(s).")]
        contigs: Vec<String>,
    },

    // Decode data written with Encode
    Decode {
        // Inputs
        // // Input fasta or fastq query file(s)
        #[arg(group = "input", required = true, help = "File with encoded fastX data.")]
        input_path: PathBuf,

        // // Prebuilt index
        #[arg(short = 'i', long = "index", help_heading = "Input", required = true, help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs")]
        index_prefix: Option<String>,
    },
    // Interactive graph visualization for Graph output
        #[cfg(feature = "graph-view")]
        GraphView {
        // Inputs
        // // Input fasta or fastq query file(s)
        #[arg(group = "input", required = true, help = "File with encoded fastX data.")]
        input_path: PathBuf,

        // // Prebuilt index
        #[arg(short = 'i', long = "index", help_heading = "Input", required = true, help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs")]
        index_prefix: Option<String>,

        #[arg(long = "color", required = true, help = "Color to plot")]
        color: usize,
    },
}

/// The input files listed in `input_list_file`: a path, or a name and a path
/// separated by `delimiter`, on each line.
pub fn read_input_list(
    input_list_file: &str,
    delimiter: u8,
) -> Result<Vec<String>, String> {
    let error = |e: &dyn std::fmt::Display| format!("reading --input-list {}: {}", input_list_file, e);
    let fs = std::fs::File::open(input_list_file).map_err(|e| error(&e))?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(false)
        .from_reader(fs);
    reader.records().map(|line| {
        let record = line.map_err(|e| error(&e))?;
        Ok(record[record.len().min(2) - 1].to_string())
    }).collect()
}

/// Parses a gap length range: `<min>-<max>`, `<min>-inf` or `none`.
fn parse_gap_range(s: &str) -> Result<ntcomp::tigs::GapRange, String> {
    if s == "none" { return Ok(None) }
    let (lo, hi) = s.split_once('-').ok_or_else(|| format!("expected <min>-<max|inf> or none, got {}", s))?;
    let lo: usize = lo.parse().map_err(|_| format!("invalid minimum in {}", s))?;
    let hi: usize = if hi == "inf" { usize::MAX } else { hi.parse().map_err(|_| format!("invalid maximum in {}", s))? };
    if lo > hi { return Err(format!("empty range {}", s)) }
    Ok(Some((lo, hi)))
}
