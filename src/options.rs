use clap::{Args as ClapArgs, Parser, Subcommand};

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
    /// Common benchmark arguments.
    #[command(flatten)]
    pub common: CommonArgs,
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
}

#[cfg(test)]
mod tests {
    use super::{Args, Command};
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
            "--ignore-abundance-error",
        ]);

        assert_eq!(args.log_level, "debug");

        match args.command {
            Command::Strain(strain_args) => {
                assert_eq!(strain_args.common.meta.as_deref(), Some("meta.tsv"));
                assert_eq!(strain_args.common.outprefix.as_deref(), Some("out"));
                assert!(strain_args.common.ignore_abundance_error);
            }
            Command::Profile(_) => panic!("Expected strain subcommand"),
        }
    }
}
