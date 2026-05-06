# geno-parse
Development Status: Early-Stage / Pre-Production. Please be aware that this software is under active development. The current implementation is an early version, and significant changes to the functionality, APIs, and overall structure should be expected in future updates.


A lightweight, memory-safe CLI tool written in Rust for parsing, filtering, and quality control of genomic data.
Designed for bioinformaticians who need fast preprocessing of sequencing data without heavyweight dependencies.

## Features

**FASTQ Processing:**
- Phred quality score calculation (ASCII 33-based)
- GC content analysis per read
- Configurable quality and length filtering
- Per-read statistics with aggregated summaries

**VCF Processing:**
- Variant classification (SNP, insertion, deletion, MNP, complex)
- QUAL-based filtering
- Per-chromosome variant distribution
- Pass/fail status tracking

**General:**
- Multi-threaded processing with Rayon
- Gzip-compressed input support
- TSV and JSON output formats
- Comprehensive error handling
- Zero unsafe code

## Installation

Requires Rust 1.56 or later. Install from source:

```bash
git clone https://github.com/ekantkannam/geno-parse.git
cd geno-parse
cargo build --release
```

The binary is available at `./target/release/geno-parse`.

## Usage

### FASTQ Quality Control

Filter reads by quality and length:

```bash
geno-parse fastq-qc -i input.fastq -q 20 -l 50 -t 8
```

With gzip input and JSON output:

```bash
geno-parse fastq-qc -i input.fastq.gz -f json -o results.json
```

Parameters:
- `-i, --input <FILE>` — Input FASTQ file
- `-q, --min-quality <N>` — Minimum mean Phred quality (default: 20)
- `-l, --min-length <N>` — Minimum read length in bp (default: 50)
- `-t, --threads <N>` — Number of threads (default: 4)
- `-f, --format <FORMAT>` — Output format: tsv or json (default: tsv)
- `-o, --output <FILE>` — Write to file (default: stdout)

### VCF Summary

Summarize and classify variants:

```bash
geno-parse vcf-summary -i variants.vcf -q 30 -f json
```

Parameters:
- `-i, --input <FILE>` — Input VCF file
- `-q, --min-qual <N>` — Minimum variant QUAL score (default: 0)
- `-f, --format <FORMAT>` — Output format: tsv or json (default: tsv)
- `-o, --output <FILE>` — Write to file (default: stdout)

### View Help

```bash
geno-parse info
geno-parse fastq-qc --help
geno-parse vcf-summary --help
```

## Example Output

### FASTQ-QC (TSV)

```
# geno-parse FASTQ-QC Summary
# Total Reads	5
# Passed Reads	3
# Failed Reads	2
# Total Bases	234
# Mean Read Length	46.80
# Mean Quality Score	32.00
# Mean GC Content	0.5089

read_id	length	mean_quality	gc_content	passed_filter
SRR001666.1	72	40.00	0.4722	PASS
SRR001666.2	72	40.00	0.4722	PASS
SRR001666.3	20	0.00	0.6000	FAIL
```

### VCF Summary (TSV)

```
# geno-parse VCF Summary
# Total Variants	10
# SNPs	6
# Insertions	3
# Deletions	1
# Mean QUAL	68.70

chrom	pos	id	ref	alt	qual	filter	variant_type
chr1	10000	rs12345	A	G	99.00	PASS	SNP
chr2	15000	rs45678	ATG	A	60.00	PASS	DEL
```

## Project Structure

```
geno-parse/
├── src/
│   ├── main.rs       — CLI entry point and command dispatch
│   ├── cli.rs        — Argument parsing with clap
│   ├── fastq.rs      — FASTQ parsing and QC logic
│   ├── vcf.rs        — VCF parsing and variant classification
│   ├── output.rs     — TSV/JSON output formatting
│   ├── errors.rs     — Error types
│   └── lib.rs        — Library exports
├── tests/
│   └── integration.rs — Integration tests
├── data/
│   ├── sample.fastq
│   └── sample.vcf
├── Cargo.toml
├── LICENSE
└── README.md
```

## Testing

Run the test suite:

```bash
cargo test
```

Run tests with output:

```bash
cargo test -- --nocapture
```

## Building for Release

Optimized binary with LTO:

```bash
cargo build --release
```

The binary is optimized for speed and minimal binary size.

## Use Cases

- **Preprocessing pipelines** — Add to Snakemake, Nextflow, or shell workflows
- **QC reporting** — Generate JSON summaries for downstream processing
- **Quick validation** — Verify FASTQ/VCF quality before expensive computations
- **Teaching** — Study genomic data processing in memory-safe Rust

## Performance

On modern hardware with 8+ cores, typical performance:
- **FASTQ**: ~1M reads/min with gzip compression
- **VCF**: ~500K variants/min
- Memory usage is proportional to input file size

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contributing

Contributions are welcome. Please ensure code passes tests and follows Rust conventions:

```bash
cargo fmt
cargo clippy
cargo test
```

## Roadmap

Potential future features:
- BAM/CRAM format support
- Variant annotation
- Statistical summaries (coverage, depth)
- HTML report generation
- Parallel file processing

