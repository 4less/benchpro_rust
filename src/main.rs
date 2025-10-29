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
use utils::time;


pub fn logo() -> String {
    let logo = "
Benchpro statistics
".to_string();
    logo
}

fn main() {
    eprintln!("{}", logo());

    let args: Args = Args::parse();
    let (duration, _) = time(|| run(&args));

    eprintln!("Benchpro took {:?}", duration);
}
