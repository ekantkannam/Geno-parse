//! VCF (Variant Call Format) parser and summarizer.
//!
//! Parses VCF files, extracting FORMAT and INFO fields, classifying variants.

use crate::errors::{GenoError, Result};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Classification of variant type
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub enum VariantType {
    Snp,
    Insertion,
    Deletion,
    Mnp,
    Complex,
    Multiallelic,
}

impl std::fmt::Display for VariantType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VariantType::Snp => write!(f, "SNP"),
            VariantType::Insertion => write!(f, "INS"),
            VariantType::Deletion => write!(f, "DEL"),
            VariantType::Mnp => write!(f, "MNP"),
            VariantType::Complex => write!(f, "COMPLEX"),
            VariantType::Multiallelic => write!(f, "MULTIALLELIC"),
        }
    }
}

/// A single VCF variant record
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub reference: String,
    pub alt: String,
    pub qual: f64,
    pub filter: String,
    pub variant_type: VariantType,
    pub info: HashMap<String, String>,
    pub gt: Option<String>,
    pub dp: Option<u32>,
    pub af: Option<f64>,
}

/// Aggregate VCF statistics
#[derive(Debug, Default, Serialize, JsonSchema)]
pub struct VcfSummary {
    pub total_variants: usize,
    pub snps: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub mnps: usize,
    pub complex: usize,
    pub multiallelic: usize,
    pub pass_filter: usize,
    pub mean_qual: f64,
    pub records: Vec<VcfRecord>,
    pub chromosomes: HashMap<String, usize>,
}

/// Determine variant type from REF and ALT alleles
fn classify_variant(reference: &str, alt: &str) -> VariantType {
    if alt.contains(',') {
        return VariantType::Multiallelic;
    }
    let ref_len = reference.len();
    let alt_len = alt.len();
    match (ref_len, alt_len) {
        (1, 1) => VariantType::Snp,
        (r, a) if r == a => VariantType::Mnp,
        (r, a) if r < a => VariantType::Insertion,
        (r, a) if r > a => VariantType::Deletion,
        _ => VariantType::Complex,
    }
}

/// Parse a QUAL field (may be '.' for missing)
fn parse_qual(s: &str) -> f64 {
    if s == "." {
        0.0
    } else {
        s.parse::<f64>().unwrap_or(0.0)
    }
}

/// Parse INFO field into HashMap
fn parse_info(info_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    if info_str == "." {
        return map;
    }
    for item in info_str.split(';') {
        let mut parts = item.splitn(2, '=');
        let key = parts.next().unwrap_or("").to_string();
        let val = parts.next().unwrap_or("").to_string();
        if !key.is_empty() {
            map.insert(key, val);
        }
    }
    map
}

/// Parse a VCF file and return a `VcfSummary`
pub fn run_vcf_summary(input: &str, min_qual: f64, benchmark: bool) -> Result<VcfSummary> {
    use std::time::Instant;
    let start_time = Instant::now();
    let mut bytes_read = 0;

    let file = File::open(input).map_err(|e| GenoError::FileIo {
        path: input.into(),
        source: e,
    })?;
    let reader = BufReader::with_capacity(64 * 1024, file);

    let mut records = Vec::new();
    let mut line_num = 0usize;

    #[allow(clippy::explicit_counter_loop)]
    for line_res in reader.lines() {
        let line = line_res.map_err(|e| GenoError::FileIo {
            path: input.into(),
            source: e,
        })?;
        bytes_read += line.len() + 1; // approx
        line_num += 1;

        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            return Err(GenoError::Parse {
                path: input.into(),
                line: line_num,
                message: format!("Expected ≥7 tab-separated fields, got {}", fields.len()),
            });
        }

        let chrom = fields[0].to_string();
        let pos = fields[1]
            .parse::<u64>()
            .map_err(|_| GenoError::InvalidVcf {
                path: input.into(),
                line: line_num,
                message: format!("Invalid POS '{}'", fields[1]),
            })?;
        let id = fields[2].to_string();
        let reference = fields[3].to_string();
        let alt = fields[4].to_string();
        let qual = parse_qual(fields[5]);
        let filter = fields[6].to_string();

        let info_str = if fields.len() > 7 { fields[7] } else { "." };
        let info = parse_info(info_str);

        let mut gt = None;
        let mut dp = None;
        let mut af = None;

        if fields.len() >= 10 {
            let format_str = fields[8];
            let sample_str = fields[9];
            let format_keys: Vec<&str> = format_str.split(':').collect();
            let sample_vals: Vec<&str> = sample_str.split(':').collect();

            for (k, v) in format_keys.iter().zip(sample_vals.iter()) {
                match *k {
                    "GT" => gt = Some(v.to_string()),
                    "DP" => dp = v.parse::<u32>().ok(),
                    "AF" => af = v.parse::<f64>().ok(),
                    _ => {}
                }
            }
        }

        if qual < min_qual {
            continue;
        }

        let variant_type = classify_variant(&reference, &alt);

        records.push(VcfRecord {
            chrom,
            pos,
            id,
            reference,
            alt,
            qual,
            filter,
            variant_type,
            info,
            gt,
            dp,
            af,
        });
    }

    let total_variants = records.len();
    let mut snps = 0;
    let mut insertions = 0;
    let mut deletions = 0;
    let mut mnps = 0;
    let mut complex = 0;
    let mut multiallelic = 0;
    let mut pass_filter = 0;
    let mut sum_qual = 0.0;

    let mut chromosomes = HashMap::new();

    for r in &records {
        match r.variant_type {
            VariantType::Snp => snps += 1,
            VariantType::Insertion => insertions += 1,
            VariantType::Deletion => deletions += 1,
            VariantType::Mnp => mnps += 1,
            VariantType::Complex => complex += 1,
            VariantType::Multiallelic => multiallelic += 1,
        }
        if r.filter == "PASS" {
            pass_filter += 1;
        }
        sum_qual += r.qual;
        *chromosomes.entry(r.chrom.clone()).or_insert(0) += 1;
    }

    let mean_qual = if total_variants == 0 {
        0.0
    } else {
        sum_qual / total_variants as f64
    };

    if benchmark {
        let elapsed = start_time.elapsed().as_secs_f64();
        let mb = bytes_read as f64 / 1_048_576.0;
        eprintln!(">> Benchmark: {:.2} seconds", elapsed);
        eprintln!(">> Speed: {:.2} MB/s", mb / elapsed);
        eprintln!(
            ">> Throughput: {:.2} variants/sec",
            total_variants as f64 / elapsed
        );
    }

    Ok(VcfSummary {
        total_variants,
        snps,
        insertions,
        deletions,
        mnps,
        complex,
        multiallelic,
        pass_filter,
        mean_qual,
        records,
        chromosomes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_multiallelic() {
        assert_eq!(classify_variant("A", "T,C"), VariantType::Multiallelic);
    }
}
