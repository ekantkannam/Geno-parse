//! Custom error types for geno-parse

use thiserror::Error;

/// All errors that can occur in the genomic parsing engine
#[derive(Debug, Error)]
pub enum GenoError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("Invalid FASTQ record: {0}")]
    InvalidFastq(String),

    #[error("Invalid VCF record: {0}")]
    InvalidVcf(String),

    #[error("Serialization error: {0}")]
    Serialization(#[from] serde_json::Error),
}

pub type Result<T> = std::result::Result<T, GenoError>;
