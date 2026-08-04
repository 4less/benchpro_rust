//! Error type for the `align` subcommand.

use std::path::PathBuf;

use thiserror::Error;

/// Anything that can go wrong while scoring an alignment benchmark.
///
/// Every variant names the file (and where possible the row or line) it came from: an alignment
/// run has hundreds of input paths, and "no such file" without one of them is not actionable.
#[derive(Debug, Error)]
pub enum AlignError {
    /// The samplesheet could not be read or has the wrong shape.
    #[error("meta file '{path}': {message}")]
    Meta {
        /// Samplesheet path.
        path: PathBuf,
        /// What is wrong with it.
        message: String,
    },

    /// The samplesheet parsed, but its rows do not describe a runnable benchmark. Carries every
    /// problem found, so one pass over the file reports all of them instead of the first.
    #[error("meta file '{path}' is not valid:\n{}", .problems.iter().map(|p| format!("  - {p}")).collect::<Vec<_>>().join("\n"))]
    MetaValidation {
        /// Samplesheet path.
        path: PathBuf,
        /// One message per problem, in row order.
        problems: Vec<String>,
    },

    /// An input file could not be read.
    #[error("cannot read '{path}': {source}")]
    Io {
        /// Offending path.
        path: PathBuf,
        /// Underlying I/O error.
        source: std::io::Error,
    },

    /// A line of a plain-text input (truth, contig2genome, `.fai`) is malformed.
    #[error("{path}:{line}: {message}")]
    Parse {
        /// Offending path.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// What was expected.
        message: String,
    },

    /// The alignment file and its truth share no read at all, so every metric would be zero.
    /// Almost always a read-id convention mismatch rather than a genuinely useless aligner.
    #[error(
        "'{alignment}' and truth '{truth}' share no read id ({aligned} primary alignments, \
         {truth_reads} truth reads). Read ids must match: the truth keys on '<id>' + mate column, \
         the alignment on QNAME with any trailing /1,/2 stripped and the mate taken from the FLAG"
    )]
    NoSharedReads {
        /// Alignment path.
        alignment: PathBuf,
        /// Truth path.
        truth: PathBuf,
        /// Primary alignments parsed.
        aligned: usize,
        /// Reads in the truth.
        truth_reads: usize,
    },

    /// Writing an output table failed.
    #[error("cannot write '{path}': {message}")]
    Output {
        /// Output path.
        path: PathBuf,
        /// Underlying error.
        message: String,
    },
}

/// Convenience alias for fallible `align` operations.
pub type AlignResult<T> = Result<T, AlignError>;

impl AlignError {
    /// Wraps an I/O error with the path it happened on.
    ///
    /// # Arguments
    ///
    /// * `path` - File the operation was attempted on
    /// * `source` - The underlying I/O error
    ///
    /// # Returns
    ///
    /// An [`AlignError::Io`] carrying both.
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }

    /// Reports a malformed line of a plain-text input.
    ///
    /// # Arguments
    ///
    /// * `path` - File being parsed
    /// * `line` - 1-based line number
    /// * `message` - What was expected on that line
    ///
    /// # Returns
    ///
    /// An [`AlignError::Parse`] naming the exact line.
    pub fn parse(path: impl Into<PathBuf>, line: usize, message: impl Into<String>) -> Self {
        Self::Parse {
            path: path.into(),
            line,
            message: message.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::AlignError;

    #[test]
    fn parse_error_names_the_line() {
        let err = AlignError::parse("truth.tsv", 3, "expected 5 tab separated fields, found 2");
        assert_eq!(
            err.to_string(),
            "truth.tsv:3: expected 5 tab separated fields, found 2"
        );
    }

    #[test]
    fn validation_error_lists_every_problem() {
        let err = AlignError::MetaValidation {
            path: "meta.tsv".into(),
            problems: vec!["row 2: no such file".into(), "row 3: duplicate".into()],
        };
        let text = err.to_string();
        assert!(text.contains("  - row 2: no such file"));
        assert!(text.contains("  - row 3: duplicate"));
    }
}
