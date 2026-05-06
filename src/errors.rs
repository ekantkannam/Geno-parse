//! Custom error types for geno-parse

use std::path::PathBuf;
use thiserror::Error;

/// All errors that can occur in the genomic parsing engine
#[derive(Debug, Error)]
pub enum GenoError {
    #[error("I/O error while reading/writing file '{path}': {source}. Check if the file exists and you have permissions.")]
    FileIo {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error in '{path}' at line {line}: {message}")]
    Parse {
        path: PathBuf,
        line: usize,
        message: String,
    },

    #[error("Invalid FASTQ format in '{path}' near line {line}: {message}. Check that the file follows standard FASTQ format.")]
    InvalidFastq {
        path: PathBuf,
        line: usize,
        message: String,
    },

    #[error("Invalid VCF format in '{path}' at line {line}: {message}. Check that the file follows standard VCF format.")]
    InvalidVcf {
        path: PathBuf,
        line: usize,
        message: String,
    },

    #[error("Serialization error: {0}")]
    Serialization(#[from] serde_json::Error),

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

pub type Result<T> = std::result::Result<T, GenoError>;
