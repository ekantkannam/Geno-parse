//! VCF (Variant Call Format) parser and summarizer

use crate::errors::{GenoError, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use serde::{Deserialize, Serialize};

/// A single VCF variant record
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub reference: String,
    pub alt: String,
    pub qual: f64,
    pub filter: String,
    pub variant_type: VariantType,
}

/// Classification of variant type
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum VariantType {
    Snp,
    Insertion,
    Deletion,
    Mnp,
    Complex,
}

impl std::fmt::Display for VariantType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VariantType::Snp => write!(f, "SNP"),
            VariantType::Insertion => write!(f, "INS"),
            VariantType::Deletion => write!(f, "DEL"),
            VariantType::Mnp => write!(f, "MNP"),
            VariantType::Complex => write!(f, "COMPLEX"),
        }
    }
}

/// Aggregate VCF statistics
#[derive(Debug, Default)]
pub struct VcfSummary {
    pub total_variants: usize,
    pub snps: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub mnps: usize,
    pub complex: usize,
    pub pass_filter: usize,
    pub mean_qual: f64,
    pub records: Vec<VcfRecord>,
    pub chromosomes: std::collections::HashMap<String, usize>,
}

/// Determine variant type from REF and ALT alleles
fn classify_variant(reference: &str, alt: &str) -> VariantType {
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

/// Parse a VCF file and return a `VcfSummary`
pub fn run_vcf_summary(input: &str, min_qual: f64) -> Result<VcfSummary> {
    let file = File::open(input)?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    let mut line_num = 0usize;

    for line in reader.lines() {
        let line = line?;
        line_num += 1;

        // Skip header and meta lines
        if line.starts_with('#') {
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 7 {
            return Err(GenoError::Parse {
                line: line_num,
                message: format!("Expected ≥7 tab-separated fields, got {}", fields.len()),
            });
        }

        let chrom = fields[0].to_string();
        let pos = fields[1].parse::<u64>().map_err(|_| GenoError::InvalidVcf(
            format!("Invalid POS '{}' at line {}", fields[1], line_num)
        ))?;
        let id = fields[2].to_string();
        let reference = fields[3].to_string();
        let alt = fields[4].to_string();
        let qual = parse_qual(fields[5]);
        let filter = fields[6].to_string();

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
        });
    }

    // Compute summary statistics
    let total_variants = records.len();
    let snps = records.iter().filter(|r| r.variant_type == VariantType::Snp).count();
    let insertions = records.iter().filter(|r| r.variant_type == VariantType::Insertion).count();
    let deletions = records.iter().filter(|r| r.variant_type == VariantType::Deletion).count();
    let mnps = records.iter().filter(|r| r.variant_type == VariantType::Mnp).count();
    let complex = records.iter().filter(|r| r.variant_type == VariantType::Complex).count();
    let pass_filter = records.iter().filter(|r| r.filter == "PASS").count();
    let mean_qual = if total_variants == 0 {
        0.0
    } else {
        records.iter().map(|r| r.qual).sum::<f64>() / total_variants as f64
    };

    let mut chromosomes = std::collections::HashMap::new();
    for rec in &records {
        *chromosomes.entry(rec.chrom.clone()).or_insert(0) += 1;
    }

    Ok(VcfSummary {
        total_variants,
        snps,
        insertions,
        deletions,
        mnps,
        complex,
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
    fn test_classify_snp() {
        assert_eq!(classify_variant("A", "T"), VariantType::Snp);
    }

    #[test]
    fn test_classify_insertion() {
        assert_eq!(classify_variant("A", "ATG"), VariantType::Insertion);
    }

    #[test]
    fn test_classify_deletion() {
        assert_eq!(classify_variant("ATG", "A"), VariantType::Deletion);
    }

    #[test]
    fn test_classify_mnp() {
        assert_eq!(classify_variant("AT", "GC"), VariantType::Mnp);
    }

    #[test]
    fn test_parse_qual_dot() {
        assert_eq!(parse_qual("."), 0.0);
    }

    #[test]
    fn test_parse_qual_number() {
        assert!((parse_qual("50.5") - 50.5).abs() < 0.001);
    }
}
