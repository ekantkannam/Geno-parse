//! FASTQ file parser with quality control and filtering

use crate::errors::{GenoError, Result};
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

/// Represents a single FASTQ read record
#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub header: String,
    pub sequence: String,
    pub quality: String,
}

/// Statistics for a single read
#[derive(Debug, Clone)]
pub struct ReadStats {
    pub header: String,
    pub length: usize,
    pub mean_quality: f64,
    pub gc_content: f64,
    pub passed_filter: bool,
}

/// Aggregate statistics over all reads
#[derive(Debug, Default)]
pub struct FastqSummary {
    pub total_reads: usize,
    pub passed_reads: usize,
    pub failed_reads: usize,
    pub total_bases: usize,
    pub mean_read_length: f64,
    pub mean_quality_score: f64,
    pub mean_gc_content: f64,
    pub read_stats: Vec<ReadStats>,
}

/// Calculate mean Phred quality score from quality string
fn mean_phred_score(quality: &str) -> f64 {
    if quality.is_empty() {
        return 0.0;
    }
    let sum: f64 = quality
        .bytes()
        .map(|b| (b.saturating_sub(33)) as f64)
        .sum();
    sum / quality.len() as f64
}

/// Calculate GC content of a sequence (0.0 to 1.0)
fn gc_content(sequence: &str) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }
    let gc = sequence
        .bytes()
        .filter(|&b| b == b'G' || b == b'C' || b == b'g' || b == b'c')
        .count();
    gc as f64 / sequence.len() as f64
}

/// Open a file, automatically handling .gz compression
fn open_reader(path: &str) -> Result<Box<dyn Read>> {
    let file = File::open(path)?;
    if path.ends_with(".gz") {
        Ok(Box::new(GzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

/// Parse FASTQ records from a reader into a Vec
fn parse_records(reader: impl Read) -> Result<Vec<FastqRecord>> {
    let buf = BufReader::new(reader);
    let mut lines = buf.lines();
    let mut records = Vec::new();
    let mut line_num = 0usize;

    loop {
        // Header line
        let header = match lines.next() {
            Some(Ok(l)) => l,
            Some(Err(e)) => return Err(GenoError::Io(e)),
            None => break,
        };
        line_num += 1;

        if header.is_empty() {
            continue;
        }
        if !header.starts_with('@') {
            return Err(GenoError::Parse {
                line: line_num,
                message: format!("Expected '@' at start of FASTQ header, got: {}", &header[..header.len().min(20)]),
            });
        }

        // Sequence
        let sequence = lines.next()
            .ok_or_else(|| GenoError::InvalidFastq("Unexpected EOF after header".into()))??;
        line_num += 1;

        // Plus separator
        let _plus = lines.next()
            .ok_or_else(|| GenoError::InvalidFastq("Missing '+' separator line".into()))??;
        line_num += 1;

        // Quality
        let quality = lines.next()
            .ok_or_else(|| GenoError::InvalidFastq("Missing quality line".into()))??;
        line_num += 1;

        if sequence.len() != quality.len() {
            return Err(GenoError::InvalidFastq(format!(
                "Sequence length {} != quality length {} at record near line {}",
                sequence.len(), quality.len(), line_num
            )));
        }

        records.push(FastqRecord {
            header: header[1..].to_string(), // strip '@'
            sequence,
            quality,
        });
    }

    Ok(records)
}

/// Run FASTQ quality control on a file
///
/// Returns a `FastqSummary` with per-read statistics and aggregate totals.
pub fn run_fastq_qc(
    input: &str,
    min_quality: f64,
    min_length: usize,
    threads: usize,
) -> Result<FastqSummary> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap_or(());

    let reader = open_reader(input)?;
    let records = parse_records(reader)?;

    if records.is_empty() {
        return Ok(FastqSummary::default());
    }

    // Process records in parallel
    let stats: Vec<ReadStats> = records
        .par_iter()
        .map(|rec| {
            let mean_quality = mean_phred_score(&rec.quality);
            let length = rec.sequence.len();
            let gc = gc_content(&rec.sequence);
            let passed = mean_quality >= min_quality && length >= min_length;
            ReadStats {
                header: rec.header.clone(),
                length,
                mean_quality,
                gc_content: gc,
                passed_filter: passed,
            }
        })
        .collect();

    let total_reads = stats.len();
    let passed_reads = stats.iter().filter(|s| s.passed_filter).count();
    let failed_reads = total_reads - passed_reads;
    let total_bases: usize = stats.iter().map(|s| s.length).sum();
    let mean_read_length = total_bases as f64 / total_reads as f64;
    let mean_quality_score = stats.iter().map(|s| s.mean_quality).sum::<f64>() / total_reads as f64;
    let mean_gc_content = stats.iter().map(|s| s.gc_content).sum::<f64>() / total_reads as f64;

    Ok(FastqSummary {
        total_reads,
        passed_reads,
        failed_reads,
        total_bases,
        mean_read_length,
        mean_quality_score,
        mean_gc_content,
        read_stats: stats,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mean_phred_score() {
        // ASCII 73 = 'I', Phred = 73-33 = 40
        let score = mean_phred_score("IIII");
        assert!((score - 40.0).abs() < 0.01);
    }

    #[test]
    fn test_gc_content() {
        assert!((gc_content("GCGC") - 1.0).abs() < 0.01);
        assert!((gc_content("ATAT") - 0.0).abs() < 0.01);
        assert!((gc_content("GCAT") - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_parse_records_valid() {
        let data = b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTT\n+\nAAAA\n";
        let records = parse_records(&data[..]).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].header, "read1");
        assert_eq!(records[0].sequence, "ACGTACGT");
    }

    #[test]
    fn test_parse_records_invalid_header() {
        let data = b"read1\nACGT\n+\nIIII\n";
        assert!(parse_records(&data[..]).is_err());
    }

    #[test]
    fn test_parse_records_length_mismatch() {
        let data = b"@read1\nACGT\n+\nII\n";
        assert!(parse_records(&data[..]).is_err());
    }
}
