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
}
