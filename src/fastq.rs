//! FASTQ file parser with quality control and filtering.
//!
//! Supports streaming architecture, paired-end processing, and adapter trimming.

use crate::errors::{GenoError, Result};
use crossbeam_channel::bounded;
use flate2::read::MultiGzDecoder;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::thread;
use std::time::Instant;

/// Represents a single FASTQ read record
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct FastqRecord {
    pub header: String,
    pub sequence: String,
    pub quality: String,
}

/// Statistics for a single read
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct ReadStats {
    pub header: String,
    pub length: usize,
    pub mean_quality: f64,
    pub gc_content: f64,
    pub passed_filter: bool,
    pub adapter_trimmed: usize,
}

/// Paired Read Statistics
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct PairedReadStats {
    pub r1: ReadStats,
    pub r2: ReadStats,
    pub both_passed: bool,
}

/// Aggregate statistics over all reads
#[derive(Debug, Default, Serialize, Deserialize, JsonSchema)]
pub struct FastqSummary {
    pub total_reads: usize,
    pub passed_reads: usize,
    pub failed_reads: usize,
    pub total_bases: usize,
    pub mean_read_length: f64,
    pub mean_quality_score: f64,
    pub mean_gc_content: f64,

    // Adapter trimming stats
    pub total_reads_trimmed: usize,
    pub total_bases_trimmed: usize,
    pub length_distribution_pre: BTreeMap<usize, usize>,
    pub length_distribution_post: BTreeMap<usize, usize>,

    // Paired stats
    pub is_paired: bool,
    pub concordance_rate: f64,
    pub both_passed: usize,
    pub one_failed: usize,

    // Individual records
    pub read_stats: Vec<ReadStats>,
    pub paired_stats: Vec<PairedReadStats>,
}

const CHUNK_SIZE: usize = 10000;
const UNIVERSAL_ADAPTER: &str = "AGATCGGAAGAGC";

/// Calculate mean Phred quality score
fn mean_phred_score(quality: &str) -> f64 {
    if quality.is_empty() {
        return 0.0;
    }
    let sum: f64 = quality.bytes().map(|b| (b.saturating_sub(33)) as f64).sum();
    sum / quality.len() as f64
}

/// Calculate GC content
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

/// Trim adapter from the 3' end using a simple suffix match
fn trim_adapter(sequence: &mut String, quality: &mut String, adapter: &str) -> usize {
    if adapter.is_empty() || sequence.len() < 3 {
        return 0;
    }
    let seq_bytes = sequence.as_bytes();
    let adapter_bytes = adapter.as_bytes();

    let min_overlap = 3;
    let max_overlap = seq_bytes.len().min(adapter_bytes.len());

    for overlap in (min_overlap..=max_overlap).rev() {
        let seq_suffix = &seq_bytes[seq_bytes.len() - overlap..];
        let adapter_prefix = &adapter_bytes[..overlap];

        if seq_suffix == adapter_prefix {
            // Match found!
            let trim_len = overlap;
            sequence.truncate(sequence.len() - trim_len);
            quality.truncate(quality.len() - trim_len);
            return trim_len;
        }
    }
    0
}

/// Open a reader, handling bgzf/gzip
fn open_reader(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|e| GenoError::FileIo {
        path: path.into(),
        source: e,
    })?;
    if path.ends_with(".gz") {
        // MultiGzDecoder handles concatenated gzip streams (like BGZF) correctly
        Ok(Box::new(BufReader::with_capacity(
            64 * 1024,
            MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(64 * 1024, file)))
    }
}

/// Parse a single FASTQ record from a lines iterator
fn parse_record(
    lines: &mut std::io::Lines<impl BufRead>,
    line_num: &mut usize,
    path: &str,
) -> Result<Option<FastqRecord>> {
    let header = match lines.next() {
        Some(Ok(l)) => l,
        Some(Err(e)) => {
            return Err(GenoError::FileIo {
                path: path.into(),
                source: e,
            })
        }
        None => return Ok(None),
    };
    *line_num += 1;

    if header.is_empty() {
        // Allow empty lines between records? Usually not standard, but we can return None or loop.
        // Let's assume strict FASTQ.
        return Err(GenoError::InvalidFastq {
            path: path.into(),
            line: *line_num,
            message: "Empty line instead of header".into(),
        });
    }
    if !header.starts_with('@') {
        return Err(GenoError::InvalidFastq {
            path: path.into(),
            line: *line_num,
            message: format!(
                "Expected '@' at start of FASTQ header, got: {}",
                &header[..header.len().min(20)]
            ),
        });
    }

    let sequence = lines
        .next()
        .ok_or_else(|| GenoError::InvalidFastq {
            path: path.into(),
            line: *line_num,
            message: "Unexpected EOF after header".into(),
        })?
        .map_err(|e| GenoError::FileIo {
            path: path.into(),
            source: e,
        })?;
    *line_num += 1;

    let _plus = lines
        .next()
        .ok_or_else(|| GenoError::InvalidFastq {
            path: path.into(),
            line: *line_num,
            message: "Missing '+' separator line".into(),
        })?
        .map_err(|e| GenoError::FileIo {
            path: path.into(),
            source: e,
        })?;
    *line_num += 1;

    let quality = lines
        .next()
        .ok_or_else(|| GenoError::InvalidFastq {
            path: path.into(),
            line: *line_num,
            message: "Missing quality line".into(),
        })?
        .map_err(|e| GenoError::FileIo {
            path: path.into(),
            source: e,
        })?;
    *line_num += 1;

    if sequence.len() != quality.len() {
        return Err(GenoError::InvalidFastq {
            path: path.into(),
            line: *line_num,
            message: format!(
                "Sequence length {} != quality length {}",
                sequence.len(),
                quality.len()
            ),
        });
    }

    Ok(Some(FastqRecord {
        header: header[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string(), // Keep only ID part of header for pairing
        sequence,
        quality,
    }))
}

enum ChunkType {
    Single(Vec<FastqRecord>),
    Paired(Vec<(FastqRecord, FastqRecord)>),
}

struct WorkerResult {
    single_stats: Vec<ReadStats>,
    paired_stats: Vec<PairedReadStats>,
    pre_lengths: BTreeMap<usize, usize>,
    post_lengths: BTreeMap<usize, usize>,
    reads_trimmed: usize,
    bases_trimmed: usize,
}

#[allow(clippy::too_many_arguments)]
pub fn run_fastq_qc(
    input: &str,
    input2: Option<&str>,
    min_quality: f64,
    min_length: usize,
    adapter1: Option<&str>,
    adapter2: Option<&str>,
    threads: usize,
    benchmark: bool,
) -> Result<FastqSummary> {
    let start_time = Instant::now();

    // For single-end
    let ad1 = adapter1.unwrap_or(UNIVERSAL_ADAPTER).to_string();
    let ad2 = adapter2.unwrap_or(UNIVERSAL_ADAPTER).to_string();

    let is_paired = input2.is_some();
    let input_path = input.to_string();
    let input2_path = input2.map(|s| s.to_string());

    let (tx_chunk, rx_chunk) = bounded::<ChunkType>(threads * 2);
    let (tx_result, rx_result) = bounded::<WorkerResult>(threads * 2);

    // Reader thread
    let reader_handle = thread::spawn(move || -> Result<usize> {
        let mut total_bytes_estimate = 0;
        if is_paired {
            let path2 = input2_path.unwrap();
            let mut lines1 = open_reader(&input_path)?.lines();
            let mut lines2 = open_reader(&path2)?.lines();
            let mut l1 = 0;
            let mut l2 = 0;

            loop {
                let mut chunk = Vec::with_capacity(CHUNK_SIZE);
                for _ in 0..CHUNK_SIZE {
                    let r1 = parse_record(&mut lines1, &mut l1, &input_path)?;
                    let r2 = parse_record(&mut lines2, &mut l2, &path2)?;

                    match (r1, r2) {
                        (Some(rec1), Some(rec2)) => {
                            // verify pairing IDs if needed, but often /1 /2 suffixes differ. We assume order is correct.
                            total_bytes_estimate += rec1.sequence.len() * 2 + 100;
                            total_bytes_estimate += rec2.sequence.len() * 2 + 100;
                            chunk.push((rec1, rec2));
                        }
                        (None, None) => break,
                        _ => {
                            return Err(GenoError::Other(anyhow::anyhow!(
                                "Paired-end files have different number of records!"
                            )))
                        }
                    }
                }
                if chunk.is_empty() {
                    break;
                }
                if tx_chunk.send(ChunkType::Paired(chunk)).is_err() {
                    break; // Receiver hung up
                }
            }
        } else {
            let mut lines = open_reader(&input_path)?.lines();
            let mut l = 0;

            loop {
                let mut chunk = Vec::with_capacity(CHUNK_SIZE);
                for _ in 0..CHUNK_SIZE {
                    if let Some(rec) = parse_record(&mut lines, &mut l, &input_path)? {
                        total_bytes_estimate += rec.sequence.len() * 2 + 100;
                        chunk.push(rec);
                    } else {
                        break;
                    }
                }
                if chunk.is_empty() {
                    break;
                }
                if tx_chunk.send(ChunkType::Single(chunk)).is_err() {
                    break;
                }
            }
        }
        Ok(total_bytes_estimate)
    });

    // Worker threads
    let mut worker_handles = Vec::new();
    for _ in 0..threads {
        let rx = rx_chunk.clone();
        let tx = tx_result.clone();
        let a1 = ad1.clone();
        let a2 = ad2.clone();

        let handle = thread::spawn(move || {
            for chunk in rx {
                let mut res = WorkerResult {
                    single_stats: Vec::new(),
                    paired_stats: Vec::new(),
                    pre_lengths: BTreeMap::new(),
                    post_lengths: BTreeMap::new(),
                    reads_trimmed: 0,
                    bases_trimmed: 0,
                };

                let mut process_rec = |mut rec: FastqRecord, adapter: &str| -> ReadStats {
                    let pre_len = rec.sequence.len();
                    *res.pre_lengths.entry(pre_len).or_insert(0) += 1;

                    let trimmed = trim_adapter(&mut rec.sequence, &mut rec.quality, adapter);
                    if trimmed > 0 {
                        res.reads_trimmed += 1;
                        res.bases_trimmed += trimmed;
                    }

                    let post_len = rec.sequence.len();
                    *res.post_lengths.entry(post_len).or_insert(0) += 1;

                    let mq = mean_phred_score(&rec.quality);
                    let gc = gc_content(&rec.sequence);
                    let passed = mq >= min_quality && post_len >= min_length;

                    ReadStats {
                        header: rec.header,
                        length: post_len,
                        mean_quality: mq,
                        gc_content: gc,
                        passed_filter: passed,
                        adapter_trimmed: trimmed,
                    }
                };

                match chunk {
                    ChunkType::Single(vec) => {
                        for rec in vec {
                            res.single_stats.push(process_rec(rec, &a1));
                        }
                    }
                    ChunkType::Paired(vec) => {
                        for (r1, r2) in vec {
                            let stat1 = process_rec(r1, &a1);
                            let stat2 = process_rec(r2, &a2);
                            let both = stat1.passed_filter && stat2.passed_filter;
                            res.paired_stats.push(PairedReadStats {
                                r1: stat1,
                                r2: stat2,
                                both_passed: both,
                            });
                        }
                    }
                }

                let _ = tx.send(res);
            }
        });
        worker_handles.push(handle);
    }

    drop(tx_result);

    let mut summary = FastqSummary {
        is_paired,
        ..Default::default()
    };

    for res in rx_result {
        summary.total_reads_trimmed += res.reads_trimmed;
        summary.total_bases_trimmed += res.bases_trimmed;

        for (k, v) in res.pre_lengths {
            *summary.length_distribution_pre.entry(k).or_insert(0) += v;
        }
        for (k, v) in res.post_lengths {
            *summary.length_distribution_post.entry(k).or_insert(0) += v;
        }

        for stat in res.single_stats {
            summary.total_reads += 1;
            summary.total_bases += stat.length;
            summary.mean_quality_score += stat.mean_quality;
            summary.mean_gc_content += stat.gc_content;
            if stat.passed_filter {
                summary.passed_reads += 1;
            } else {
                summary.failed_reads += 1;
            }
            summary.read_stats.push(stat);
        }

        for pstat in res.paired_stats {
            // Count R1 and R2 as separate reads for overall totals
            summary.total_reads += 2;
            summary.total_bases += pstat.r1.length + pstat.r2.length;
            summary.mean_quality_score += pstat.r1.mean_quality + pstat.r2.mean_quality;
            summary.mean_gc_content += pstat.r1.gc_content + pstat.r2.gc_content;

            if pstat.both_passed {
                summary.both_passed += 1;
                summary.passed_reads += 2;
            } else {
                summary.one_failed += 1;
                if pstat.r1.passed_filter {
                    summary.passed_reads += 1;
                } else {
                    summary.failed_reads += 1;
                }
                if pstat.r2.passed_filter {
                    summary.passed_reads += 1;
                } else {
                    summary.failed_reads += 1;
                }
            }

            summary.paired_stats.push(pstat);
        }
    }

    // Wait for reader and get stats
    let bytes_estimate = reader_handle.join().unwrap()?;
    for handle in worker_handles {
        handle.join().unwrap();
    }

    if summary.total_reads > 0 {
        summary.mean_read_length = summary.total_bases as f64 / summary.total_reads as f64;
        summary.mean_quality_score /= summary.total_reads as f64;
        summary.mean_gc_content /= summary.total_reads as f64;
    }

    let total_pairs = summary.both_passed + summary.one_failed;
    if total_pairs > 0 {
        summary.concordance_rate = summary.both_passed as f64 / total_pairs as f64;
    }

    if benchmark {
        let elapsed = start_time.elapsed().as_secs_f64();
        let mb = bytes_estimate as f64 / 1_048_576.0;
        eprintln!(">> Benchmark: {:.2} seconds", elapsed);
        eprintln!(">> Speed: {:.2} MB/s", mb / elapsed);
        eprintln!(
            ">> Throughput: {:.2} reads/sec",
            summary.total_reads as f64 / elapsed
        );
    }

    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_adapter() {
        let mut seq = "ACGTACGTAGATCGGAAGAGC".to_string();
        let mut qual = "IIIIIIIIIIIIIIIIIIIII".to_string();
        let trimmed = trim_adapter(&mut seq, &mut qual, UNIVERSAL_ADAPTER);
        assert_eq!(trimmed, 13);
        assert_eq!(seq, "ACGTACGT");
        assert_eq!(qual, "IIIIIIII");

        // Partial match at end
        let mut seq = "ACGTACGTAGAT".to_string();
        let mut qual = "IIIIIIIIIIII".to_string();
        let trimmed = trim_adapter(&mut seq, &mut qual, UNIVERSAL_ADAPTER);
        assert_eq!(trimmed, 4); // "AGAT" matches
        assert_eq!(seq, "ACGTACGT");
    }
}
