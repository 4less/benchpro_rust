//! Read-alignment benchmark.
//!
//! Where [`crate::benchpro`] asks *which taxa are in a sample*, this module asks a different
//! question about a different input: given a read aligner's SAM/PAF output and per-read ground
//! truth, *did each read go to the right place, and is the alignment it reports any good?*
//!
//! It shares none of the taxonomic machinery (no `Taxonomy`, no `TaxonomicRank`) — only benchpro's
//! shape: a meta samplesheet in, aggregated TSVs out.
//!
//! The module scores files that already exist; it never runs an aligner and never measures time.
//! Orchestration (building indexes, timing runs, cache control) belongs to the harness that calls
//! it.

pub mod cigar;
pub mod error;
pub mod meta;
pub mod sam;
pub mod truth;

use log::info;

use crate::options::AlignArgs;
use error::AlignResult;
use meta::AlignMeta;

/// Runs the alignment benchmark described by `args.meta`.
///
/// # Arguments
///
/// * `args` - Parsed `align` subcommand arguments
///
/// # Returns
///
/// `Ok(())` once every output table has been written.
///
/// # Errors
///
/// Returns an [`error::AlignError`] when the samplesheet is invalid, an input cannot be read or an
/// output cannot be written.
pub fn run_align(args: &AlignArgs) -> AlignResult<()> {
    let meta = AlignMeta::from_path(&args.meta, args)?;
    info!(
        "Loaded {} row(s) from '{}': {} dataset(s), {} tool(s)",
        meta.rows.len(),
        args.meta.display(),
        meta.datasets().len(),
        meta.tools().len()
    );

    if args.validate_meta {
        info!("Meta file is valid");
        return Ok(());
    }

    Ok(())
}
