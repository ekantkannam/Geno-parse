.PHONY: build test clean benchmark bench-prep

build:
	cargo build --release

test:
	cargo test

clean:
	cargo clean
	rm -f data/synthetic_1M.fastq data/synthetic_100K.vcf tools/generate_synthetic

tools/generate_synthetic: tools/generate_synthetic.rs
	rustc -O tools/generate_synthetic.rs -o tools/generate_synthetic

bench-prep: tools/generate_synthetic build
	@mkdir -p data
	@echo "Generating synthetic 1M-read FASTQ..."
	@./tools/generate_synthetic fastq data/synthetic_1M.fastq
	@echo "Generating synthetic 100K-variant VCF..."
	@./tools/generate_synthetic vcf data/synthetic_100K.vcf

benchmark: bench-prep
	@echo "========================================"
	@echo "           BENCHMARK RESULTS            "
	@echo "========================================"
	@echo "Hardware: $$(uname -sm)"
	@echo "Threads: 8"
	@echo ""
	
	@echo "[ geno-parse (FASTQ) ]"
	@if command -v /usr/bin/time > /dev/null; then \
		/usr/bin/time -l ./target/release/geno-parse fastq-qc -i data/synthetic_1M.fastq -t 8 --benchmark > /dev/null; \
	else \
		time ./target/release/geno-parse fastq-qc -i data/synthetic_1M.fastq -t 8 --benchmark > /dev/null; \
	fi
	@echo ""

	@echo "[ fastp ]"
	@if command -v fastp > /dev/null; then \
		if command -v /usr/bin/time > /dev/null; then \
			/usr/bin/time -l fastp -i data/synthetic_1M.fastq -w 8 -o /dev/null > /dev/null 2>&1; \
		else \
			time fastp -i data/synthetic_1M.fastq -w 8 -o /dev/null > /dev/null 2>&1; \
		fi \
	else \
		echo "fastp not installed, skipping..."; \
	fi
	@echo ""

	@echo "[ seqkit stats (FASTQ) ]"
	@if command -v seqkit > /dev/null; then \
		if command -v /usr/bin/time > /dev/null; then \
			/usr/bin/time -l seqkit stats -j 8 data/synthetic_1M.fastq > /dev/null 2>&1; \
		else \
			time seqkit stats -j 8 data/synthetic_1M.fastq > /dev/null 2>&1; \
		fi \
	else \
		echo "seqkit not installed, skipping..."; \
	fi
	@echo ""

	@echo "[ geno-parse (VCF) ]"
	@if command -v /usr/bin/time > /dev/null; then \
		/usr/bin/time -l ./target/release/geno-parse vcf-summary -i data/synthetic_100K.vcf --benchmark > /dev/null; \
	else \
		time ./target/release/geno-parse vcf-summary -i data/synthetic_100K.vcf --benchmark > /dev/null; \
	fi
	@echo ""
