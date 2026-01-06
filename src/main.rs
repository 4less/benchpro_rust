#![feature(adt_const_params)]
#![feature(generic_arg_infer)]
#![feature(let_chains)]
#![feature(iter_collect_into)]
#[macro_use]
pub mod common;
pub mod profile;
pub mod format;
pub mod test_data;
pub mod options;
pub mod utils;
pub mod benchpro;
pub mod meta;
pub mod taxonomy;
pub mod profile_handler;
pub mod tree_handler;
pub mod tree_adjusted_benchmarks;
pub mod taxonomic_profiling;

use benchpro::run;
use clap::Parser;
use options::Args;
use env_logger::Env;
use log::info;
use utils::time;


/// Return the Benchpro banner string.
///
/// # Returns
///
/// Static banner shown at startup.
pub fn logo() -> String {
    let logo = "
Benchpro statistics
".to_string();
    logo
}

fn main() {
    let args: Args = Args::parse();
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

    info!("{}", logo());
    let (duration, _) = time(|| run(&args));

    info!("Benchpro took {:?}", duration);
}
