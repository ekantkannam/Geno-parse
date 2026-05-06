//! Integration tests for geno-parse

use std::io::Write;
use tempfile::NamedTempFile;

fn write_temp_fastq(content: &str) -> NamedTempFile {
    let mut f = NamedTempFile::new().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

fn write_temp_vcf(content: &str) -> NamedTempFile {
    let mut f = NamedTempFile::new().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f
}

#[test]
fn test_fastq_qc_all_pass() {
    // Reads with quality 'I' = Phred 40, length 8 — should all pass
    let data = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nGCATGCAT\n+\nIIIIIIII\n";
    let f = write_temp_fastq(data);
    let summary = geno_parse::fastq::run_fastq_qc(f.path().to_str().unwrap(), 20.0, 4, 2).unwrap();
    assert_eq!(summary.total_reads, 2);
    assert_eq!(summary.passed_reads, 2);
    assert_eq!(summary.failed_reads, 0);
}

#[test]
fn test_fastq_qc_filter_by_quality() {
    // First read has quality '!' = Phred 0, should fail min_quality=20
    let data = "@r1\nACGT\n+\n!!!!\n@r2\nACGTACGT\n+\nIIIIIIII\n";
    let f = write_temp_fastq(data);
    let summary = geno_parse::fastq::run_fastq_qc(f.path().to_str().unwrap(), 20.0, 4, 2).unwrap();
    assert_eq!(summary.total_reads, 2);
    assert_eq!(summary.passed_reads, 1);
    assert_eq!(summary.failed_reads, 1);
}

#[test]
fn test_fastq_qc_filter_by_length() {
    let data = "@r1\nAA\n+\nII\n@r2\nACGTACGT\n+\nIIIIIIII\n";
    let f = write_temp_fastq(data);
    let summary = geno_parse::fastq::run_fastq_qc(f.path().to_str().unwrap(), 20.0, 5, 2).unwrap();
    assert_eq!(summary.passed_reads, 1);
}

#[test]
fn test_vcf_summary_basic() {
    let vcf_data = "##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        chr1\t100\t.\tA\tT\t50\tPASS\t.\n\
        chr1\t200\t.\tATG\tA\t30\tPASS\t.\n\
        chr2\t300\t.\tC\tCTT\t40\t.\t.\n";
    let f = write_temp_vcf(vcf_data);
    let summary = geno_parse::vcf::run_vcf_summary(f.path().to_str().unwrap(), 0.0).unwrap();
    assert_eq!(summary.total_variants, 3);
    assert_eq!(summary.snps, 1);
    assert_eq!(summary.deletions, 1);
    assert_eq!(summary.insertions, 1);
    assert_eq!(summary.pass_filter, 2);
}

#[test]
fn test_vcf_min_qual_filter() {
    let vcf_data = "##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        chr1\t100\t.\tA\tT\t10\tPASS\t.\n\
        chr1\t200\t.\tG\tC\t60\tPASS\t.\n";
    let f = write_temp_vcf(vcf_data);
    let summary = geno_parse::vcf::run_vcf_summary(f.path().to_str().unwrap(), 30.0).unwrap();
    assert_eq!(summary.total_variants, 1);
    assert_eq!(summary.snps, 1);
}
