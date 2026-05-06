//! geno-parse: Fast genomic data parsing and quality control
//!
//! A lightweight CLI tool for FASTQ and VCF processing in genomic analysis workflows.
//! Routes subcommands to appropriate processing modules and manages output formatting.

mod cli;
mod errors;
mod fastq;
mod output;
mod vcf;

use clap::Parser;
use cli::{Cli, Commands};
use std::fs::File;
use std::io::{self, BufWriter, Write};

/// Creates an output writer that writes to either a file or stdout.
fn get_writer(output: &Option<String>) -> Box<dyn Write> {
    match output {
        Some(path) => {
            let file = File::create(path).expect("Failed to create output file");
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(io::stdout())),
    }
}

/// Displays the application banner and basic usage information.
fn print_banner() {
    println!(
        r#"
  ╔═══════════════════════════════════════════════════╗
  ║         geno-parse v0.1.0  ⚡ Genomics CLI         ║
  ║        FASTQ and VCF Parsing & Quality Control   ║
  ╚═══════════════════════════════════════════════════╝
"#
    );
}

fn main() {
    // Suppress SIGPIPE panics (e.g., when piped to `head`)
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    let cli = Cli::parse();

    if cli.schema {
        output::print_schemas();
        std::process::exit(0);
    }

    if let Some(cmd) = cli.command {
        match cmd {
            Commands::FastqQc {
                input,
                input2,
                output,
                threads,
                min_quality,
                min_length,
                adapter1,
                adapter2,
                benchmark,
                format,
            } => {
                eprintln!("⚡ Running FASTQ-QC on: {}", input);
                if let Some(ref i2) = input2 {
                    eprintln!("   Paired with: {}", i2);
                }
                eprintln!(
                    "   Threads: {}  |  Min Quality: {}  |  Min Length: {}",
                    threads, min_quality, min_length
                );

                match fastq::run_fastq_qc(
                    &input,
                    input2.as_deref(),
                    min_quality,
                    min_length,
                    adapter1.as_deref(),
                    adapter2.as_deref(),
                    threads,
                    benchmark,
                ) {
                    Ok(summary) => {
                        eprintln!(
                            "✅ Done: {}/{} reads passed ({} failed)",
                            summary.passed_reads, summary.total_reads, summary.failed_reads
                        );
                        let mut writer = get_writer(&output);
                        if let Err(e) = output::write_fastq_output(&summary, &format, &mut writer) {
                            if let errors::GenoError::Io(ref source) = e {
                                if source.kind() == std::io::ErrorKind::BrokenPipe {
                                    std::process::exit(0);
                                }
                            }
                            eprintln!("❌ Output error: {}", e);
                            std::process::exit(1);
                        }
                        if let Err(e) = writer.flush() {
                            if e.kind() == std::io::ErrorKind::BrokenPipe {
                                std::process::exit(0);
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("❌ Error: {}", e);
                        std::process::exit(1);
                    }
                }
            }

            Commands::VcfSummary {
                input,
                output,
                format,
                min_qual,
                benchmark,
            } => {
                eprintln!("⚡ Running VCF-Summary on: {}", input);
                eprintln!("   Min QUAL: {}", min_qual);

                match vcf::run_vcf_summary(&input, min_qual, benchmark) {
                    Ok(summary) => {
                        eprintln!(
                            "✅ Done: {} variants found (SNPs: {}, INS: {}, DEL: {})",
                            summary.total_variants,
                            summary.snps,
                            summary.insertions,
                            summary.deletions
                        );
                        let mut writer = get_writer(&output);
                        if let Err(e) = output::write_vcf_output(&summary, &format, &mut writer) {
                            if let errors::GenoError::Io(ref source) = e {
                                if source.kind() == std::io::ErrorKind::BrokenPipe {
                                    std::process::exit(0);
                                }
                            }
                            eprintln!("❌ Output error: {}", e);
                            std::process::exit(1);
                        }
                        if let Err(e) = writer.flush() {
                            if e.kind() == std::io::ErrorKind::BrokenPipe {
                                std::process::exit(0);
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("❌ Error: {}", e);
                        std::process::exit(1);
                    }
                }
            }

            Commands::Info => {
                print_banner();
                println!("  Subcommands:");
                println!("    fastq-qc      Quality control & filtering of FASTQ files");
                println!("    vcf-summary   Parse and summarize VCF variant files");
                println!();
                println!("  Features:");
                println!("    ✓ Multi-threaded parallel processing (Rayon & Crossbeam)");
                println!("    ✓ Gzip-compressed input (.gz) support via BGZF");
                println!("    ✓ TSV and JSON output formats with robust schema");
                println!("    ✓ Phred quality scoring & Adapter Trimming");
                println!("    ✓ GC content calculation & Paired-End metrics");
                println!("    ✓ SNP / INDEL / MNP variant classification");
                println!();
                println!("  Usage examples:");
                println!("    geno-parse fastq-qc -i sample_R1.fastq -I sample_R2.fastq -q 20 -l 50 -t 8");
                println!("    geno-parse fastq-qc -i sample.fastq.gz -f json -o results.json");
                println!("    geno-parse vcf-summary -i variants.vcf -q 30 -f tsv");
            }
        }
    } else {
        print_banner();
        println!("Run `geno-parse --help` for usage instructions.");
    }
}
