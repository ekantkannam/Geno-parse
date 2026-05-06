//! Output formatting — TSV and JSON writers

use crate::cli::OutputFormat;
use crate::errors::Result;
use crate::fastq::FastqSummary;
use crate::vcf::VcfSummary;
use schemars::schema_for;
use std::io::Write;

/// Write FASTQ QC results to a writer in the chosen format
pub fn write_fastq_output(
    summary: &FastqSummary,
    format: &OutputFormat,
    writer: &mut dyn Write,
) -> Result<()> {
    match format {
        OutputFormat::Tsv => {
            writeln!(writer, "# geno-parse FASTQ-QC Summary")?;
            writeln!(writer, "# Is Paired-End\t{}", summary.is_paired)?;
            writeln!(writer, "# Total Reads\t{}", summary.total_reads)?;
            writeln!(writer, "# Passed Reads\t{}", summary.passed_reads)?;
            writeln!(writer, "# Failed Reads\t{}", summary.failed_reads)?;
            writeln!(writer, "# Total Bases\t{}", summary.total_bases)?;
            writeln!(
                writer,
                "# Total Reads Trimmed\t{}",
                summary.total_reads_trimmed
            )?;
            writeln!(
                writer,
                "# Total Bases Trimmed\t{}",
                summary.total_bases_trimmed
            )?;
            if summary.is_paired {
                writeln!(
                    writer,
                    "# Concordance Rate\t{:.4}",
                    summary.concordance_rate
                )?;
                writeln!(writer, "# Both Passed\t{}", summary.both_passed)?;
                writeln!(writer, "# One Failed\t{}", summary.one_failed)?;
            }
            writeln!(
                writer,
                "# Mean Read Length\t{:.2}",
                summary.mean_read_length
            )?;
            writeln!(
                writer,
                "# Mean Quality Score\t{:.2}",
                summary.mean_quality_score
            )?;
            writeln!(writer, "# Mean GC Content\t{:.4}", summary.mean_gc_content)?;
            writeln!(writer)?;
            writeln!(
                writer,
                "read_id\tlength\tmean_quality\tgc_content\tpassed_filter\tadapter_trimmed"
            )?;
            for stat in &summary.read_stats {
                writeln!(
                    writer,
                    "{}\t{}\t{:.2}\t{:.4}\t{}\t{}",
                    stat.header,
                    stat.length,
                    stat.mean_quality,
                    stat.gc_content,
                    if stat.passed_filter { "PASS" } else { "FAIL" },
                    stat.adapter_trimmed
                )?;
            }
            for pstat in &summary.paired_stats {
                writeln!(
                    writer,
                    "{}/1\t{}\t{:.2}\t{:.4}\t{}\t{}",
                    pstat.r1.header,
                    pstat.r1.length,
                    pstat.r1.mean_quality,
                    pstat.r1.gc_content,
                    if pstat.r1.passed_filter {
                        "PASS"
                    } else {
                        "FAIL"
                    },
                    pstat.r1.adapter_trimmed
                )?;
                writeln!(
                    writer,
                    "{}/2\t{}\t{:.2}\t{:.4}\t{}\t{}",
                    pstat.r2.header,
                    pstat.r2.length,
                    pstat.r2.mean_quality,
                    pstat.r2.gc_content,
                    if pstat.r2.passed_filter {
                        "PASS"
                    } else {
                        "FAIL"
                    },
                    pstat.r2.adapter_trimmed
                )?;
            }
        }
        OutputFormat::Json => {
            writeln!(writer, "{}", serde_json::to_string_pretty(&summary)?)?;
        }
    }
    Ok(())
}

/// Write VCF summary results to a writer in the chosen format
pub fn write_vcf_output(
    summary: &VcfSummary,
    format: &OutputFormat,
    writer: &mut dyn Write,
) -> Result<()> {
    match format {
        OutputFormat::Tsv => {
            writeln!(writer, "# geno-parse VCF Summary")?;
            writeln!(writer, "# Total Variants\t{}", summary.total_variants)?;
            writeln!(writer, "# SNPs\t{}", summary.snps)?;
            writeln!(writer, "# Insertions\t{}", summary.insertions)?;
            writeln!(writer, "# Deletions\t{}", summary.deletions)?;
            writeln!(writer, "# MNPs\t{}", summary.mnps)?;
            writeln!(writer, "# Complex\t{}", summary.complex)?;
            writeln!(writer, "# Multiallelic\t{}", summary.multiallelic)?;
            writeln!(writer, "# PASS Filter\t{}", summary.pass_filter)?;
            writeln!(writer, "# Mean QUAL\t{:.2}", summary.mean_qual)?;
            writeln!(writer)?;
            writeln!(
                writer,
                "chrom\tpos\tid\tref\talt\tqual\tfilter\tvariant_type\tgt\tdp\taf"
            )?;
            for rec in &summary.records {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}",
                    rec.chrom,
                    rec.pos,
                    rec.id,
                    rec.reference,
                    rec.alt,
                    rec.qual,
                    rec.filter,
                    rec.variant_type,
                    rec.gt.as_deref().unwrap_or("."),
                    rec.dp.map_or(".".to_string(), |v| v.to_string()),
                    rec.af.map_or(".".to_string(), |v| format!("{:.4}", v))
                )?;
            }
        }
        OutputFormat::Json => {
            writeln!(writer, "{}", serde_json::to_string_pretty(&summary)?)?;
        }
    }
    Ok(())
}

/// Print the JSON schemas for the output structures to stdout
pub fn print_schemas() {
    let fastq_schema = schema_for!(FastqSummary);
    let vcf_schema = schema_for!(VcfSummary);

    println!("--- FASTQ Summary Schema ---");
    println!("{}", serde_json::to_string_pretty(&fastq_schema).unwrap());
    println!();
    println!("--- VCF Summary Schema ---");
    println!("{}", serde_json::to_string_pretty(&vcf_schema).unwrap());
}
