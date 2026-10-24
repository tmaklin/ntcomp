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

use clap::{Args, Parser, Subcommand};

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
    Build(BuildArgs),
    // Choose the index k-mers from a collection and write them as tigs
    Tigs(TigsArgs),
    // Encode fastX data using an SBWT index
    Encode(EncodeArgs),
    // Encode a collection of fastX files using an SBWT index
    Graph(GraphArgs),
    // Decode data written with Encode
    Decode(DecodeArgs),
    #[cfg(feature = "graph-view")]
    // Interactive graph visualization for Graph output
    GraphView(GraphViewArgs),
}

// Build the compression dictionary
#[derive(Args)]
pub struct BuildArgs {
    // Input fasta or fastq sequence file(s)
    #[arg(group = "input", required = true, help = "Sequence data file(s).")]
    pub seq_files: Vec<String>,

    #[arg(short = 'l', long = "input-list", group = "input", required = true, help_heading = "Input", help = "File with paths or tab separated name and path on each line.")]
    pub input_list: Option<String>,

    // Outputs
    #[arg(short = 'o', long = "output-prefix", required = true, help_heading = "Output", help = "Prefix for output files <prefix>.sbwt and <prefix>.lcs.")]
    pub output_prefix: String,

    // Build parameters
    // // k-mer size
    #[arg(short = 'k', default_value_t = 31, help_heading = "Build options", help = "k-mer size, larger values are slower and use more space.")]
    pub kmer_size: usize,
    // // prefix precalc
    #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options", help = "Length of precalculated prefixes included in the index.")]
    pub prefix_precalc: usize,
    // // deduplicate k-mer batches
    #[arg(short = 'd', long = "dedup-batches", default_value_t = false, help_heading = "Build options", help = "Deduplicate k-mer batches to save some memory.")]
    pub dedup_batches: bool,

    // Resources
    // // Threads
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    pub num_threads: usize,
    // // Memory in GB
    #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options", help = "Memory available when building on temp disk space (in gigabytes).")]
    pub mem_gb: usize,
    // // Temporary directory
    #[arg(long = "temp-dir", required = false, help_heading = "Build options", help = "Build on temporary disk space at this path instead of in-memory.")]
    pub temp_dir: Option<String>,

    // Verbosity
    #[arg(long = "verbose", default_value_t = false)]
    pub verbose: bool,
}

// Choose the index k-mers from a collection and write them as tigs
#[derive(Args)]
pub struct TigsArgs {
    // Input fasta or fastq sequence file(s)
    #[arg(group = "input", required = true, help = "Assembly file(s), one per input; the tigs are seeded in this order.")]
    pub seq_files: Vec<String>,

    #[arg(short = 'l', long = "input-list", group = "input", required = true, help_heading = "Input", help = "File with paths or tab separated name and path on each line.")]
    pub input_list: Option<String>,

    // Outputs
    #[arg(short = 'o', long = "output-prefix", required = true, help_heading = "Output", help = "Prefix for <prefix>.tigs.tac (and the outputs below).")]
    pub output_prefix: String,

    #[arg(long = "fasta", default_value_t = false, help_heading = "Output", help = "Also write the tigs to <prefix>.tigs.fasta.")]
    pub fasta: bool,

    #[arg(long = "index", default_value_t = false, help_heading = "Output", help = "Also write the index built from the tigs to <prefix>.sbwt and <prefix>.lcs.")]
    pub index: bool,

    // Tig parameters
    #[arg(short = 'k', default_value_t = 255, help_heading = "Tig options", help = "k-mer size.")]
    pub kmer_size: usize,

    #[arg(short = 'T', long = "min-assemblies", default_value_t = 5, help_heading = "Tig options", help = "Core: canonical k-mers present in at least this many inputs.")]
    pub min_assemblies: u32,

    #[arg(long = "gap-join", default_value = "1-inf", value_parser = parse_gap_range, help_heading = "Tig options", help = "Lengths (in k-mers) of gaps between two core tig ends to fill: <min>-<max|inf> or none.")]
    pub gap_join: ntcomp::tigs::GapRange,

    #[arg(long = "gap-one-sided", default_value = "1-inf", value_parser = parse_gap_range, help_heading = "Tig options", help = "Lengths of gaps with exactly one flank a core tig end to fill.")]
    pub gap_one_sided: ntcomp::tigs::GapRange,

    #[arg(long = "gap-bubble", default_value = "350-inf", value_parser = parse_gap_range, help_heading = "Tig options", help = "Lengths of gaps with neither flank a core tig end to fill.")]
    pub gap_bubble: ntcomp::tigs::GapRange,

    #[arg(long = "gap-tip", default_value = "250-inf", value_parser = parse_gap_range, help_heading = "Tig options", help = "Lengths of runs at contig ends to fill.")]
    pub gap_tip: ntcomp::tigs::GapRange,

    #[arg(long = "gap-keep-repairable", default_value_t = false, help_heading = "Tig options", help = "Also fill the bubbles that SNP and indel editing repairs against the core (left out by default).")]
    pub gap_keep_repairable: bool,

    #[arg(long = "orient-sampling", default_value_t = 8, help_heading = "Tig options", help = "Sample one in this many 31-mers for orienting the tigs.")]
    pub orient_sampling: u64,

    #[arg(long = "partitions", default_value_t = 16, help_heading = "Tig options", help = "Count the core in this many hash partitions (more use less memory but read the inputs more times).")]
    pub partitions: u64,

    // Index build parameters
    #[arg(short = 'p', long = "prefix-precalc", default_value_t = 8, help_heading = "Build options", help = "Length of precalculated prefixes included in the index.")]
    pub prefix_precalc: usize,

    #[arg(short = 'm', long = "mem-gb", default_value_t = 4, help_heading = "Build options", help = "Memory for building the indexes on temporary disk space (in gigabytes).")]
    pub mem_gb: usize,

    #[arg(long = "temp-dir", required = false, help_heading = "Build options", help = "Temporary disk space for building the indexes [default: <prefix>.tmp].")]
    pub temp_dir: Option<String>,

    // Resources
    #[arg(short = 't', long = "threads", default_value_t = 1, help = "Number of threads.")]
    pub num_threads: usize,

    #[arg(long = "verbose", default_value_t = false)]
    pub verbose: bool,
}

// Encode fastX data using an SBWT index
#[derive(Args)]
pub struct EncodeArgs {
    // Inputs
    // // Input fasta or fastq query file(s)
    #[arg(group = "input", required = true, help = "Query file with sequence data.")]
    pub query_file: PathBuf,

    // // Prebuilt index
    #[arg(short = 'i', long = "index", help_heading = "Input", required = true, help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs")]
    pub index_prefix: String,

}

// Encode a collection of fastX files using an SBWT index
#[derive(Args)]
pub struct GraphArgs {
    // Inputs
    // // Input fasta or fastq query file(s)
    #[arg(group = "input", required = true, help = "Querys file with sequence data.")]
    pub query_files: Vec<PathBuf>,

    // // Prebuilt index
    #[arg(short = 'i', long = "index", help_heading = "Input", required = true, help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs")]
    pub index_prefix: String,

    #[arg(long = "fasta", default_value_t = false, help = ".fasta output")]
    pub fasta: bool,

    #[arg(long = "fasta-cols", default_value_t = 80, help = ".fasta columns")]
    pub fasta_columns: usize,

    // // Collapse linear chains
    #[arg(long = "collapse-chains", default_value_t = false, help = "Store linear chains of nodes as single edges.")]
    pub collapse_chains: bool,

    // // Choose contig strands
    #[arg(long = "reorient", default_value_t = false, help = "Store every contig on the strand that shares the most with the contigs before it.")]
    pub reorient: bool,

    #[arg(long = "reorient-sampling", default_value_t = 32, value_parser = clap::value_parser!(u64).range(1..), help = "Sample one in this many 31-mers for choosing the strands.")]
    pub reorient_sampling: u64,

    // // Edit SNPs towards the index
    #[arg(long = "edit-snps", default_value_t = false, help = "Edit SNPs and small indels against the index before encoding (reads the input twice).")]
    pub edit_snps: bool,


    // Resources
    // // Threads
    #[arg(short = 't', long = "threads", default_value_t = 1, help = "Number of threads.")]
    pub num_threads: usize,

    // // Batch size
    #[arg(long = "batch-size", default_value_t = 64, value_parser = clap::value_parser!(u64).range(1..), help = "Read and process this many input files at a time (when decoding: this many files' worth of contigs).")]
    pub batch_size: u64,

    // Decompression
    // // Toggle decompression
    #[arg(short = 'd', long = "decompress", default_value_t = false, help = "Force decompression", group = "mode")]
    pub decompress: bool,
    // // Retrieve specific contig(s) when decompressing
    #[arg(long = "color", help_heading = "Decompression", help = "Decode by contig index(es) (0-based).")]
    pub colors: Vec<u32>,

    #[arg(long = "contig", help_heading = "Decompression", help = "Decode by contig name(s).")]
    pub contigs: Vec<String>,
}

// Decode data written with Encode
#[derive(Args)]
pub struct DecodeArgs {
    // Inputs
    // // Input fasta or fastq query file(s)
    #[arg(group = "input", required = true, help = "File with encoded fastX data.")]
    pub input_path: PathBuf,

    // // Prebuilt index
    #[arg(short = 'i', long = "index", help_heading = "Input", required = true, help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs")]
    pub index_prefix: String,
}

// Interactive graph visualization for Graph output
#[cfg(feature = "graph-view")]
#[derive(Args)]
pub struct GraphViewArgs {
    // Inputs
    // // Input fasta or fastq query file(s)
    #[arg(group = "input", required = true, help = "File with encoded fastX data.")]
    pub input_path: PathBuf,

    // // Prebuilt index
    #[arg(short = 'i', long = "index", help_heading = "Input", required = true, help = "Prefix for prebuilt <prefix>.sbwt and <prefix>.lcs")]
    pub index_prefix: String,

    #[arg(long = "color", required = true, help = "Color to plot")]
    pub color: usize,
}

/// The input files: `files`, then those listed in `input_list`.
fn input_files(
    files: &[String],
    input_list: &Option<String>,
) -> Result<Vec<String>, String> {
    let mut in_files = files.to_vec();
    if let Some(list) = input_list {
        in_files.extend(read_input_list(list, b'\t')?);
    }
    Ok(in_files)
}

impl BuildArgs {
    /// The input files: the arguments, then those in `--input-list`.
    pub fn input_files(&self) -> Result<Vec<String>, String> {
        input_files(&self.seq_files, &self.input_list)
    }

    /// The options for building the index.
    pub fn build_opts(&self) -> ntcomp::index::BuildOpts {
        ntcomp::index::BuildOpts {
            k: self.kmer_size,
            num_threads: self.num_threads,
            prefix_precalc: self.prefix_precalc,
            mem_gb: self.mem_gb,
            dedup_batches: self.dedup_batches,
            temp_dir: self.temp_dir.clone(),
        }
    }
}

impl TigsArgs {
    /// The input files: the arguments, then those in `--input-list`.
    pub fn input_files(&self) -> Result<Vec<String>, String> {
        input_files(&self.seq_files, &self.input_list)
    }

    /// The directory for building the indexes on disk (in memory, sorting
    /// the k-mers takes k / 4 bytes each): `--temp-dir`, or
    /// `<prefix>.tmp`, which is to be removed afterwards (true).
    pub fn temp_dir(&self) -> (String, bool) {
        match &self.temp_dir {
            Some(dir) => (dir.clone(), false),
            None => (format!("{}.tmp", self.output_prefix), true),
        }
    }

    /// The options for choosing the tigs, building the indexes in `temp_dir`.
    pub fn tig_opts(
        &self,
        temp_dir: &str,
    ) -> ntcomp::tigs::TigOpts {
        ntcomp::tigs::TigOpts {
            min_assemblies: self.min_assemblies,
            partitions: self.partitions,
            gaps: ntcomp::tigs::GapRules {
                join: self.gap_join,
                one_sided: self.gap_one_sided,
                bubble: self.gap_bubble,
                tip: self.gap_tip,
                drop_repairable_bubbles: !self.gap_keep_repairable,
            },
            orient_sampling: self.orient_sampling,
            build_opts: ntcomp::index::BuildOpts {
                k: self.kmer_size,
                num_threads: self.num_threads,
                prefix_precalc: self.prefix_precalc,
                mem_gb: self.mem_gb,
                temp_dir: Some(temp_dir.to_string()),
                ..Default::default()
            },
        }
    }
}

impl GraphArgs {
    /// The options for writing the archive.
    pub fn archive_opts(&self) -> ntcomp::archive::ArchiveOpts {
        ntcomp::archive::ArchiveOpts {
            collapse_chains: self.collapse_chains,
            reorient: self.reorient,
            reorient_sampling: self.reorient_sampling,
            edit: self.edit_snps,
        }
    }

    /// Contigs decoded at a time: about `--batch-size` input files' worth of
    /// an archive of `n_contigs` contigs from `n_files` files.
    pub fn chunk_size(
        &self,
        n_contigs: u32,
        n_files: usize,
    ) -> usize {
        self.batch_size as usize * (n_contigs as usize).div_ceil(n_files.max(1)).max(1)
    }
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
