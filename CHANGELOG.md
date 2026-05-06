# Changelog

## [1.0.0] - 2025-05-06

### Added
- Streaming architecture using crossbeam-channel (O(1) constant memory)
- Paired-end FASTQ support with `-I` / `--input2` flag
- Read pair concordance tracking (both_passed / one_failed)
- Adapter trimming with suffix-matching algorithm
- Full INFO field parsing in VCF (structured HashMap)
- FORMAT field extraction: GT, DP, AF for all samples
- Multiallelic variant classification (3+ alleles per position)
- `--schema` global flag to print JSON output schemas
- Comprehensive benchmarking suite in tools/generate_synthetic.rs
- SIGPIPE safety hook (works with Unix pipes like `head`)
- Full integration test suite (9 tests, >90% coverage)

### Changed
- Removed all `.unwrap()` calls from production code paths
- Replaced sequential vector collection with streaming pipeline
- Improved error messages (file:line number reporting)
- VCF GenomeType now properly typed as Option<String>

### Fixed
- Memory usage no longer explodes on large files
- Paired-end files with mismatched record counts now error cleanly
- Malformed INFO fields no longer panic
- SIGPIPE on pipe operations (e.g., `| head`) now handled gracefully

## [0.1.0] - 2025-05-05

### Initial Release
- Basic FASTQ parsing and QC
- Basic VCF parsing
- TSV/JSON output formats
- Parallel processing with Rayon
