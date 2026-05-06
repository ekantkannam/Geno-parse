//! CLI definitions using clap derive macros

use clap::{Parser, Subcommand, ValueEnum};

/// geno-parse: High-Performance Genomic Parsing Engine
#[derive(Parser, Debug)]
#[command(
    name = "geno-parse",
    version = "0.1.0",
    author = "Genomic Engine",
    about = "⚡ Fast, memory-safe CLI for parsing & QC of genomic data files",
    long_about = "A high-performance genomic parsing engine for FASTQ, VCF, and BAM files.\nBuilt with Rust for maximum throughput and minimal memory usage."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

/// Available subcommands
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Quality control and filtering of FASTQ files
    FastqQc {
        /// Input FASTQ file (plain or .gz compressed)
        #[arg(short, long, value_name = "FILE")]
        input: String,

        /// Output file (default: stdout)
        #[arg(short, long, value_name = "FILE")]
        output: Option<String>,

        /// Number of threads to use
        #[arg(short, long, default_value = "4", value_name = "N")]
        threads: usize,

        /// Minimum mean Phred quality score to keep a read
        #[arg(short = 'q', long, default_value = "20", value_name = "N")]
        min_quality: f64,

        /// Minimum read length to keep
        #[arg(short = 'l', long, default_value = "50", value_name = "N")]
        min_length: usize,

        /// Output format
        #[arg(short, long, default_value = "tsv", value_name = "FORMAT")]
        format: OutputFormat,
    },

    /// Parse and summarize VCF variant files
    VcfSummary {
        /// Input VCF file
        #[arg(short, long, value_name = "FILE")]
        input: String,

        /// Output file (default: stdout)
        #[arg(short, long, value_name = "FILE")]
        output: Option<String>,

        /// Output format
        #[arg(short, long, default_value = "tsv", value_name = "FORMAT")]
        format: OutputFormat,

        /// Minimum variant QUAL score to include
        #[arg(short = 'q', long, default_value = "0", value_name = "N")]
        min_qual: f64,
    },

    /// Display project info and benchmark stats
    Info,
}

/// Output format options
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum OutputFormat {
    Tsv,
    Json,
}
