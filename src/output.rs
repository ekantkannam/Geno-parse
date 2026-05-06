//! Output formatting — TSV and JSON writers

use crate::cli::OutputFormat;
use crate::errors::Result;
use crate::fastq::FastqSummary;
use crate::vcf::VcfSummary;
use serde_json::json;
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
            writeln!(writer, "# Total Reads\t{}", summary.total_reads)?;
            writeln!(writer, "# Passed Reads\t{}", summary.passed_reads)?;
            writeln!(writer, "# Failed Reads\t{}", summary.failed_reads)?;
            writeln!(writer, "# Total Bases\t{}", summary.total_bases)?;
            writeln!(writer, "# Mean Read Length\t{:.2}", summary.mean_read_length)?;
            writeln!(writer, "# Mean Quality Score\t{:.2}", summary.mean_quality_score)?;
            writeln!(writer, "# Mean GC Content\t{:.4}", summary.mean_gc_content)?;
            writeln!(writer)?;
            writeln!(writer, "read_id\tlength\tmean_quality\tgc_content\tpassed_filter")?;
            for stat in &summary.read_stats {
                writeln!(
                    writer,
                    "{}\t{}\t{:.2}\t{:.4}\t{}",
                    stat.header,
                    stat.length,
                    stat.mean_quality,
                    stat.gc_content,
                    if stat.passed_filter { "PASS" } else { "FAIL" }
                )?;
            }
        }
        OutputFormat::Json => {
            let read_stats_json: Vec<_> = summary.read_stats.iter().map(|s| {
                json!({
                    "read_id": s.header,
                    "length": s.length,
                    "mean_quality": format!("{:.2}", s.mean_quality).parse::<f64>().unwrap_or(s.mean_quality),
                    "gc_content": format!("{:.4}", s.gc_content).parse::<f64>().unwrap_or(s.gc_content),
                    "passed_filter": s.passed_filter,
                })
            }).collect();

            let out = json!({
                "summary": {
                    "total_reads": summary.total_reads,
                    "passed_reads": summary.passed_reads,
                    "failed_reads": summary.failed_reads,
                    "total_bases": summary.total_bases,
                    "mean_read_length": format!("{:.2}", summary.mean_read_length).parse::<f64>().unwrap_or(summary.mean_read_length),
                    "mean_quality_score": format!("{:.2}", summary.mean_quality_score).parse::<f64>().unwrap_or(summary.mean_quality_score),
                    "mean_gc_content": format!("{:.4}", summary.mean_gc_content).parse::<f64>().unwrap_or(summary.mean_gc_content),
                },
                "reads": read_stats_json,
            });
            writeln!(writer, "{}", serde_json::to_string_pretty(&out)?)?;
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
            writeln!(writer, "# PASS Filter\t{}", summary.pass_filter)?;
            writeln!(writer, "# Mean QUAL\t{:.2}", summary.mean_qual)?;
            writeln!(writer)?;
            writeln!(writer, "chrom\tpos\tid\tref\talt\tqual\tfilter\tvariant_type")?;
            for rec in &summary.records {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
                    rec.chrom, rec.pos, rec.id, rec.reference, rec.alt, rec.qual, rec.filter, rec.variant_type
                )?;
            }
        }
        OutputFormat::Json => {
            let records_json: Vec<_> = summary.records.iter().map(|r| {
                json!({
                    "chrom": r.chrom,
                    "pos": r.pos,
                    "id": r.id,
                    "ref": r.reference,
                    "alt": r.alt,
                    "qual": r.qual,
                    "filter": r.filter,
                    "variant_type": r.variant_type.to_string(),
                })
            }).collect();

            let mut chrom_counts: Vec<_> = summary.chromosomes.iter().collect();
            chrom_counts.sort_by_key(|(k, _)| k.as_str());

            let out = json!({
                "summary": {
                    "total_variants": summary.total_variants,
                    "snps": summary.snps,
                    "insertions": summary.insertions,
                    "deletions": summary.deletions,
                    "mnps": summary.mnps,
                    "complex": summary.complex,
                    "pass_filter": summary.pass_filter,
                    "mean_qual": format!("{:.2}", summary.mean_qual).parse::<f64>().unwrap_or(summary.mean_qual),
                    "chromosomes": summary.chromosomes,
                },
                "variants": records_json,
            });
            writeln!(writer, "{}", serde_json::to_string_pretty(&out)?)?;
        }
    }
    Ok(())
}
