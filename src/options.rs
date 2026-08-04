use clap::{Args as ClapArgs, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

/// ASCII logo shown at startup and in help output.
pub const LOGO: &str = "
┓      ┓       
┣┓┏┓┏┓┏┣┓┏┓┏┓┏┓
┗┛┗ ┛┗┗┛┗┣┛┛ ┗┛
         ┛     
";

/// CLI arguments for Benchpro.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(arg_required_else_help(true))]
#[command(max_term_width = 120)]
pub struct Args {
    /// Log level: debug, info, warning, error
    #[arg(long, default_value = "info", global = true)]
    pub log_level: String,

    /// Subcommand to run.
    #[command(subcommand)]
    pub command: Command,
}

/// Supported Benchpro subcommands.
#[derive(Subcommand, Debug)]
pub enum Command {
    /// Benchmark at the profile level (current behavior).
    Profile(ProfileArgs),
    /// Benchmark at the strain level.
    Strain(StrainArgs),
    /// Merge profiles into an abundance matrix.
    Merge(MergeArgs),
    /// Compute MSA summary statistics.
    Msa(MsaArgs),
    /// Normalize bacterial profiles to a standard format.
    Normalize(NormalizeArgs),
}

/// CLI arguments for the `profile` subcommand.
#[derive(ClapArgs, Debug, Clone)]
pub struct ProfileArgs {
    /// Common benchmark arguments.
    #[command(flatten)]
    pub common: CommonArgs,
}

/// CLI arguments for the `strain` subcommand.
#[derive(ClapArgs, Debug, Clone)]
pub struct StrainArgs {
    /// Samplesheet meta file with columns: ID, Species, Meta (required); Tree, MSA, Partition, GoldMSA, GoldTree, Tool (optional).
    #[arg(short = 'm', long, required = true, value_name = "PATH")]
    pub meta: std::path::PathBuf,

    /// Output prefix for strain benchmarking files.
    /// E.g. `--outprefix results/run1` produces `results/run1.monophyly.tsv`.
    #[arg(short = 'o', long, required = true, value_name = "PREFIX")]
    pub outprefix: String,

    /// Exclude samples whose coverage (from the species meta file) is below this threshold.
    #[arg(long, value_name = "FLOAT")]
    pub cov_filter: Option<f64>,

    /// Midpoint-root each tree before computing monophyly scores.
    /// Use this when input trees are unrooted (trifurcating root), which is typical
    /// for IQ-TREE / RAxML output, to avoid LCA artefacts from the arbitrary root.
    #[arg(long, default_value_t = false)]
    pub midpoint_root: bool,
}

/// CLI arguments for the `normalize` subcommand.
#[derive(ClapArgs, Debug, Clone)]
pub struct NormalizeArgs {
    /// Input profile file to normalize
    #[arg(short = 'i', long)]
    pub input: String,

    /// Output file for normalized profile
    #[arg(short = 'o', long, required_unless_present = "detect")]
    pub output: Option<String>,

    /// Only detect and print profile format/tool/version.
    #[arg(long, default_value_t = false)]
    pub detect: bool,

    /// Output format for normalized content.
    #[arg(long, value_enum, default_value_t = NormalizeOutputFormat::Plain)]
    pub output_format: NormalizeOutputFormat,
}

/// Output format for the `normalize` subcommand.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum NormalizeOutputFormat {
    Plain,
    Cami,
}

/// Shared benchmark arguments used by multiple subcommands.
#[derive(ClapArgs, Debug, Clone)]
pub struct CommonArgs {
    /// Meta file
    #[arg(short = 'm', long)]
    pub meta: Option<String>,

    /// Output prefix for all files. E.g. '--outprefix output_folder/output' will produce
    /// the files output_folder/output_summary_stats.tsv and output_folder/output_sample_stats.tsv.
    /// Optional when using --validate-meta.
    #[arg(short = 'o', long)]
    pub outprefix: Option<String>,

    /// Enable adjusted benchmarks (tree-based adjustments)
    #[arg(long, default_value_t = false)]
    pub adjusted: bool,

    /// Allow alternative taxon matching
    #[arg(long, default_value_t = false)]
    pub allow_alternatives: bool,

    /// Validate the meta file (columns and file paths) and exit
    #[arg(long, default_value_t = false)]
    pub validate_meta: bool,

    /// Ignore abundance normalization errors (log warnings instead of crashing)
    #[arg(long, default_value_t = false)]
    pub ignore_abundance_error: bool,

    /// Ignore CAMI lineage length mismatches (use last lineage token as name)
    #[arg(long, default_value_t = false)]
    pub cami_ignore_lineage_error: bool,
}

/// CLI arguments for the `merge` subcommand.
#[derive(ClapArgs, Debug, Clone)]
pub struct MergeArgs {
    /// Input profile files (supports shell globs).
    #[arg(long, num_args = 1.., required = true, value_name = "PROFILE")]
    pub input: Vec<PathBuf>,

    /// Output TSV file.
    #[arg(long, required = true, value_name = "PATH")]
    pub output: PathBuf,

    /// Target rank to include (e.g., species, genus).
    #[arg(long, default_value = "species", value_name = "RANK")]
    pub rank: String,

    /// Ignore CAMI lineage length mismatches (use last lineage token as name)
    #[arg(long, default_value_t = false)]
    pub cami_ignore_lineage_error: bool,

    /// Regex to derive sample names from input paths using capture groups.
    #[arg(long, num_args = 1.., value_name = "REGEX")]
    pub sample_regex: Vec<String>,

    /// Print inferred sample name and path for each input and exit.
    #[arg(long, default_value_t = false)]
    pub test_sample_regex: bool,
}

/// CLI arguments for the `msa` subcommand.
#[derive(ClapArgs, Debug, Clone)]
pub struct MsaArgs {
    /// Meta file with columns: ID, MSA, Meta.
    #[arg(long, required = true, value_name = "PATH")]
    pub meta: PathBuf,

    /// Output prefix for MSA statistics files.
    #[arg(long, required = true, value_name = "PATH")]
    pub output: PathBuf,
}

#[cfg(test)]
mod tests {
    use super::{Args, Command, NormalizeOutputFormat};
    use clap::Parser;

    #[test]
    fn parses_profile_subcommand() {
        let args = Args::parse_from([
            "benchpro",
            "profile",
            "--meta",
            "meta.tsv",
            "--outprefix",
            "out",
        ]);

        match args.command {
            Command::Profile(profile_args) => {
                assert_eq!(profile_args.common.meta.as_deref(), Some("meta.tsv"));
                assert_eq!(profile_args.common.outprefix.as_deref(), Some("out"));
            }
            Command::Strain(_) => panic!("Expected profile subcommand"),
            Command::Merge(_) => panic!("Expected profile subcommand"),
            Command::Msa(_) => panic!("Expected profile subcommand"),
            Command::Normalize(_) => panic!("Expected profile subcommand"),
        }
    }

    #[test]
    fn parses_strain_subcommand() {
        let args = Args::parse_from([
            "benchpro",
            "strain",
            "--meta",
            "meta.tsv",
            "--outprefix",
            "out",
            "--log-level",
            "debug",
        ]);

        assert_eq!(args.log_level, "debug");

        match args.command {
            Command::Strain(strain_args) => {
                assert_eq!(strain_args.meta.to_string_lossy(), "meta.tsv");
                assert_eq!(strain_args.outprefix, "out");
            }
            Command::Profile(_) => panic!("Expected strain subcommand"),
            Command::Merge(_) => panic!("Expected strain subcommand"),
            Command::Msa(_) => panic!("Expected strain subcommand"),
            Command::Normalize(_) => panic!("Expected strain subcommand"),
        }
    }

    #[test]
    fn parses_merge_subcommand() {
        let args = Args::parse_from([
            "benchpro",
            "merge",
            "--input",
            "a.profile",
            "b.profile",
            "--output",
            "abundance.tsv",
            "--rank",
            "genus",
            "--cami-ignore-lineage-error",
            "--sample-regex",
            "example1",
            "example2",
        ]);

        match args.command {
            Command::Merge(merge_args) => {
                assert_eq!(merge_args.input.len(), 2);
                assert_eq!(merge_args.output.to_string_lossy(), "abundance.tsv");
                assert_eq!(merge_args.rank, "genus");
                assert!(merge_args.cami_ignore_lineage_error);
                assert_eq!(
                    merge_args.sample_regex,
                    vec!["example1".to_string(), "example2".to_string()]
                );
            }
            _ => panic!("Expected merge subcommand"),
        }
    }

    #[test]
    fn parses_msa_subcommand() {
        let args = Args::parse_from([
            "benchpro",
            "msa",
            "--meta",
            "msa_meta.tsv",
            "--output",
            "stats.tsv",
        ]);

        match args.command {
            Command::Msa(msa_args) => {
                assert_eq!(msa_args.meta.to_string_lossy(), "msa_meta.tsv");
                assert_eq!(msa_args.output.to_string_lossy(), "stats.tsv");
            }
            _ => panic!("Expected msa subcommand"),
        }
    }

    #[test]
    fn parses_normalize_subcommand() {
        let args = Args::parse_from([
            "benchpro",
            "normalize",
            "--input",
            "input.profile",
            "--output",
            "output.tsv",
        ]);

        match args.command {
            Command::Normalize(normalize_args) => {
                assert_eq!(normalize_args.input, "input.profile");
                assert_eq!(normalize_args.output.as_deref(), Some("output.tsv"));
                assert!(!normalize_args.detect);
                assert_eq!(normalize_args.output_format, NormalizeOutputFormat::Plain);
            }
            _ => panic!("Expected normalize subcommand"),
        }
    }

    #[test]
    fn parses_normalize_detect_only_subcommand() {
        let args = Args::parse_from([
            "benchpro",
            "normalize",
            "--input",
            "input.profile",
            "--detect",
        ]);

        match args.command {
            Command::Normalize(normalize_args) => {
                assert_eq!(normalize_args.input, "input.profile");
                assert_eq!(normalize_args.output, None);
                assert!(normalize_args.detect);
                assert_eq!(normalize_args.output_format, NormalizeOutputFormat::Plain);
            }
            _ => panic!("Expected normalize subcommand"),
        }
    }

    #[test]
    fn parses_normalize_cami_output_format() {
        let args = Args::parse_from([
            "benchpro",
            "normalize",
            "--input",
            "input.profile",
            "--output",
            "output.profile",
            "--output-format",
            "cami",
        ]);

        match args.command {
            Command::Normalize(normalize_args) => {
                assert_eq!(normalize_args.output_format, NormalizeOutputFormat::Cami);
            }
            _ => panic!("Expected normalize subcommand"),
        }
    }
}
