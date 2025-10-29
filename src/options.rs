use clap::Parser;


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
// #[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[command(max_term_width = 120)] // term_width sets it fixed, max term_width can be smaller
pub struct Args {
    /// Meta file
    #[arg(short = 'm', long)] // String::default()
    pub meta: Option<String>,

    /// Output prefix for all files. E.g. '--outprefix output_folder/output' will produce the files output_folder/output_summary_stats.tsv and output_folder/output_sample_stats.tsv 
    #[arg(short = 'm', long)] // String::default()
    pub outprefix: String,

    // /// Reverse read of pair (.fastq, .fq)
    // #[arg(num_args(0..), short = '2', long, default_values_t = ["".to_string()], action = clap::ArgAction::Append)]
    // pub rev: Vec<String>,
}