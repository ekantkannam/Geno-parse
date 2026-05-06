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
    let data = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nGCATGCAT\n+\nIIIIIIII\n";
    let f = write_temp_fastq(data);
    let summary = geno_parse::fastq::run_fastq_qc(
        f.path().to_str().unwrap(),
        None,
        20.0,
        4,
        None,
        None,
        2,
        false,
    )
    .unwrap();
    assert_eq!(summary.total_reads, 2);
    assert_eq!(summary.passed_reads, 2);
    assert_eq!(summary.failed_reads, 0);
}

#[test]
fn test_fastq_qc_filter_by_quality() {
    let data = "@r1\nACGT\n+\n!!!!\n@r2\nACGTACGT\n+\nIIIIIIII\n";
    let f = write_temp_fastq(data);
    let summary = geno_parse::fastq::run_fastq_qc(
        f.path().to_str().unwrap(),
        None,
        20.0,
        4,
        None,
        None,
        2,
        false,
    )
    .unwrap();
    assert_eq!(summary.total_reads, 2);
    assert_eq!(summary.passed_reads, 1);
    assert_eq!(summary.failed_reads, 1);
}

#[test]
fn test_fastq_qc_filter_by_length() {
    let data = "@r1\nAA\n+\nII\n@r2\nACGTACGT\n+\nIIIIIIII\n";
    let f = write_temp_fastq(data);
    let summary = geno_parse::fastq::run_fastq_qc(
        f.path().to_str().unwrap(),
        None,
        20.0,
        5,
        None,
        None,
        2,
        false,
    )
    .unwrap();
    assert_eq!(summary.passed_reads, 1);
}

#[test]
fn test_fastq_qc_paired_end() {
    let r1 = "@read1/1\nACGTACGT\n+\nIIIIIIII\n";
    let r2 = "@read1/2\nACGTACGT\n+\nIIIIIIII\n";
    let f1 = write_temp_fastq(r1);
    let f2 = write_temp_fastq(r2);
    let summary = geno_parse::fastq::run_fastq_qc(
        f1.path().to_str().unwrap(),
        Some(f2.path().to_str().unwrap()),
        20.0,
        4,
        None,
        None,
        2,
        false,
    )
    .unwrap();
    assert_eq!(summary.total_reads, 2);
    assert_eq!(summary.both_passed, 1);
    assert!(summary.is_paired);
    assert_eq!(summary.concordance_rate, 1.0);
}

#[test]
fn test_fastq_qc_adapter_trimming() {
    let data = "@r1\nACGTACGTAGAT\n+\nIIIIIIIIIIII\n";
    let f = write_temp_fastq(data);
    let summary = geno_parse::fastq::run_fastq_qc(
        f.path().to_str().unwrap(),
        None,
        20.0,
        4,
        None,
        None,
        2,
        false,
    )
    .unwrap();
    assert_eq!(summary.total_reads_trimmed, 1);
    assert_eq!(summary.total_bases_trimmed, 4); // "AGAT" should be trimmed
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
    let summary = geno_parse::vcf::run_vcf_summary(f.path().to_str().unwrap(), 0.0, false).unwrap();
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
    let summary =
        geno_parse::vcf::run_vcf_summary(f.path().to_str().unwrap(), 30.0, false).unwrap();
    assert_eq!(summary.total_variants, 1);
    assert_eq!(summary.snps, 1);
}

#[test]
fn test_vcf_multiallelic() {
    let vcf_data = "##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
        chr1\t100\t.\tA\tT,C\t50\tPASS\tDP=100;AF=0.5,0.5\tGT:DP\t1/2:100\n";
    let f = write_temp_vcf(vcf_data);
    let summary = geno_parse::vcf::run_vcf_summary(f.path().to_str().unwrap(), 0.0, false).unwrap();
    assert_eq!(summary.total_variants, 1);
    assert_eq!(summary.multiallelic, 1);

    let rec = &summary.records[0];
    assert_eq!(rec.variant_type, geno_parse::vcf::VariantType::Multiallelic);
    assert_eq!(rec.gt, Some("1/2".to_string()));
    assert_eq!(rec.dp, Some(100));
}

#[test]
fn test_vcf_missing_qual() {
    let vcf_data = "##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        chr1\t100\t.\tA\tT\t.\tPASS\t.\n";
    let f = write_temp_vcf(vcf_data);
    let summary = geno_parse::vcf::run_vcf_summary(f.path().to_str().unwrap(), 0.0, false).unwrap();
    assert_eq!(summary.total_variants, 1);
    assert_eq!(summary.records[0].qual, 0.0);
}
