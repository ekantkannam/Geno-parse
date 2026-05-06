//! geno-parse — High-Performance Genomic Parsing Engine
//! Entry point: routes CLI subcommands to appropriate modules.

mod cli;
mod errors;
mod fastq;
mod output;
mod vcf;

use clap::Parser;
use cli::{Cli, Commands};
use std::fs::File;
use std::io::{self, BufWriter, Write};

fn get_writer(output: &Option<String>) -> Box<dyn Write> {
    match output {
        Some(path) => {
            let file = File::create(path).expect("Failed to create output file");
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(io::stdout())),
    }
}

fn print_banner() {
    println!(r#"
  ╔═══════════════════════════════════════════════════╗
  ║        geno-parse v0.1.0  ⚡ Genomic Engine       ║
  ║   High-Performance FASTQ · VCF · BAM Processing  ║
  ╚═══════════════════════════════════════════════════╝
"#);
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::FastqQc {
            input,
            output,
            threads,
            min_quality,
            min_length,
            format,
        } => {
            eprintln!("⚡ Running FASTQ-QC on: {}", input);
            eprintln!("   Threads: {}  |  Min Quality: {}  |  Min Length: {}", threads, min_quality, min_length);

            match fastq::run_fastq_qc(&input, min_quality, min_length, threads) {
                Ok(summary) => {
                    eprintln!(
                        "✅ Done: {}/{} reads passed ({} failed)",
                        summary.passed_reads, summary.total_reads, summary.failed_reads
                    );
                    let mut writer = get_writer(&output);
                    if let Err(e) = output::write_fastq_output(&summary, &format, &mut writer) {
                        eprintln!("❌ Output error: {}", e);
                        std::process::exit(1);
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
        } => {
            eprintln!("⚡ Running VCF-Summary on: {}", input);
            eprintln!("   Min QUAL: {}", min_qual);

            match vcf::run_vcf_summary(&input, min_qual) {
                Ok(summary) => {
                    eprintln!(
                        "✅ Done: {} variants found (SNPs: {}, INS: {}, DEL: {})",
                        summary.total_variants, summary.snps, summary.insertions, summary.deletions
                    );
                    let mut writer = get_writer(&output);
                    if let Err(e) = output::write_vcf_output(&summary, &format, &mut writer) {
                        eprintln!("❌ Output error: {}", e);
                        std::process::exit(1);
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
            println!("    ✓ Multi-threaded parallel processing (Rayon)");
            println!("    ✓ Gzip-compressed input (.gz) support");
            println!("    ✓ TSV and JSON output formats");
            println!("    ✓ Phred quality scoring");
            println!("    ✓ GC content calculation");
            println!("    ✓ SNP / INDEL / MNP variant classification");
            println!();
            println!("  Usage examples:");
            println!("    geno-parse fastq-qc -i sample.fastq -q 20 -l 50 -t 8");
            println!("    geno-parse fastq-qc -i sample.fastq.gz -f json -o results.json");
            println!("    geno-parse vcf-summary -i variants.vcf -q 30 -f tsv");
        }
    }
}
