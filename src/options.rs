use clap::Parser;


/// CLI arguments for Benchpro.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(arg_required_else_help(true))]
#[command(max_term_width = 120)]
pub struct Args {
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

    /// Log level: debug, info, warning, error
    #[arg(long, default_value = "info")]
    pub log_level: String,

}
