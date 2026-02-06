#![feature(adt_const_params)]
#![feature(generic_arg_infer)]
#![feature(let_chains)]
#![feature(iter_collect_into)]
#[macro_use]
pub mod common;
pub mod benchpro;
pub mod format;
pub mod meta;
pub mod options;
pub mod profile;
pub mod profile_handler;
pub mod taxonomic_profiling;
pub mod taxonomy;
pub mod test_data;
pub mod tree_adjusted_benchmarks;
pub mod tree_handler;
pub mod utils;

use benchpro::run;
use clap::{CommandFactory, FromArgMatches, Parser};
use env_logger::Env;
use log::info;
use options::{Args, Command, LOGO};
use utils::time;

fn main() {
    let mut command = Args::command();
    command = command.before_help(LOGO);
    for subcommand in command.get_subcommands_mut() {
        *subcommand = subcommand.clone().before_help(LOGO);
    }
    let matches = command.get_matches();
    let args = Args::from_arg_matches(&matches).expect("Argument parsing failed");
    let level = match args.log_level.to_lowercase().as_str() {
        "debug" => "debug",
        "info" => "info",
        "warn" | "warning" => "warn",
        "error" => "error",
        other => {
            panic!(
                "Invalid --log-level '{}'. Use: debug, info, warning, error.",
                other
            );
        }
    };
    let env = Env::default().default_filter_or(level);
    env_logger::Builder::from_env(env).init();

    info!("{}", LOGO);
    let (duration, _) = time(|| match &args.command {
        Command::Profile(profile_args) => {
            run(&profile_args.common);
        }
        Command::Strain(strain_args) => {
            run(&strain_args.common);
        }
    });

    info!("Benchpro took {:?}", duration);
}
