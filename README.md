# ⚡ geno-parse — High-Performance Genomic Parsing Engine

`geno-parse` is a lightweight, Rust-based CLI tool for rapid parsing, filtering, and quality control of genomic data.
It is designed for life-science workflows that need fast FASTQ and VCF preprocessing without a heavyweight bioinformatics stack.

## 🚀 What this tool does

- Performs **FASTQ quality control** with Phred scoring, GC content calculation, and length-based filtering
- Parses and summarizes **VCF variant files** with classification into SNP, insertion, deletion, MNP, and complex events
- Supports **TSV and JSON** outputs for integration with analysis pipelines
- Reads **gzip-compressed FASTQ** files natively
- Uses **multi-threaded processing** via Rayon for faster throughput

## ✅ Features

- FASTQ-QC:
  - per-read mean Phred quality score
  - GC content calculation
  - configurable minimum quality and length thresholds
- VCF summary:
  - QUAL-based filtering
  - variant type classification
  - pass/fail counts and per-chromosome statistics
- Output formats:
  - `tsv`
  - `json`
- Minimal dependency surface with safe Rust code
- Tested with unit and integration tests

## 📦 Installation

Build the project with Cargo:

```bash
cargo build --release
```

The binary will be available at `./target/release/geno-parse`.

## 🧪 Run tests

```bash
cargo test
```

## 💻 Usage

### FASTQ Quality Control

```bash
./target/release/geno-parse fastq-qc \
  -i sample.fastq \
  -q 20 \
  -l 50 \
  -t 8
```

```bash
./target/release/geno-parse fastq-qc \
  -i sample.fastq.gz \
  -f json \
  -o results.json
```

### VCF Summary

```bash
./target/release/geno-parse vcf-summary \
  -i variants.vcf \
  -q 30
```

```bash
./target/release/geno-parse vcf-summary \
  -i variants.vcf \
  -f json \
  -o variants_summary.json
```

### Info

```bash
./target/release/geno-parse info
```

## 📁 Sample output

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
SRR001666.1_length=72	72	40.00	0.4722	PASS
SRR001666.3_LOW_QUALITY_length=20	20	0.00	0.6000	FAIL
```

### VCF Summary (TSV)

```
# geno-parse VCF Summary
# Total Variants	10
# SNPs	6
# Insertions	3
# Deletions	1
# MNPs	0
# Complex	0
# PASS Filter	8
# Mean QUAL	68.70

chrom	pos	id	ref	alt	qual	filter	variant_type
chr1	10000	rs12345	A	G	99.00	PASS	SNP
chr2	15000	rs45678	ATG	A	60.00	PASS	DEL
```

## 🧩 Project structure

```text
geno-parse/
├── Cargo.toml
├── README.md
├── data/
│   ├── sample.fastq
│   └── sample.vcf
├── src/
│   ├── cli.rs
│   ├── errors.rs
│   ├── fastq.rs
│   ├── lib.rs
│   ├── main.rs
│   ├── output.rs
│   └── vcf.rs
└── tests/
    └── integration.rs
```

## 🧠 Architecture

- `src/main.rs` — CLI entry point and command dispatch
- `src/cli.rs` — clap-based argument parsing
- `src/fastq.rs` — FASTQ parsing, QC, and read metrics
- `src/vcf.rs` — VCF parsing and variant summarization
- `src/output.rs` — TSV/JSON formatting and writer logic
- `src/errors.rs` — custom error definitions
- `tests/integration.rs` — functional integration tests

## 📌 Publishing to GitHub

1. Create a repository on GitHub named `geno-parse`.
2. Initialize git in the project folder if needed:
   ```bash
git init
git add .
git commit -m "Initial geno-parse implementation"
git branch -M main
git remote add origin https://github.com/<your-user>/geno-parse.git
git push -u origin main
```
3. Add a `LICENSE` file if you want to publish under an open-source license.
4. Add release notes or GitHub Actions workflow later if desired.

## ✨ Notes

- This tool is intended as a **preprocessing utility**, not a full bioinformatics pipeline.
- Use it to add FASTQ/QC and VCF summary steps to larger workflows like Snakemake, Nextflow, or bash scripts.

## 📣 Want to improve it?

- Add support for **BAM/CRAM** input
- Add **variant annotation** support
- Add **read-level filtering output** to a separate file
- Add **summary dashboard** or HTML report generation
