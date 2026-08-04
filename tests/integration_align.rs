//! End-to-end tests for `benchpro align`.
//!
//! The fixture in `data/test_data/align/` is deliberately tiny and hand-checkable — every expected
//! number below is derived from how it was constructed, not recorded from a run:
//!
//! Reference: three markers, `1001_geneA` and `1001_geneB` on `genomeA`, `2002_geneA` on
//! `genomeB`, 120 bp each. Truth: 10 pairs (20 mates), all really from `1001_geneA`.
//!
//! | tool | what it does |
//! |---|---|
//! | `good` (SAM) | pairs 1-8 exact at the true locus; pair 9 on the *other* marker of the right species; pair 10 not placed at all |
//! | `sloppy` (SAM) | pairs 1-5 right but with no `NM` tag; pairs 6-8 on the wrong genome; pair 9 mate 1 malformed and mate 2 preceded by a secondary record; pair 10 mate 1 unmapped, mate 2 soft clipped 5 bp off |
//! | `mapper` (PAF) | pairs 1-9 at the true locus, keeping the `/1`,`/2` suffix SAM tools strip |

use std::collections::HashMap;
use std::{fs, path::PathBuf, process::Command};

mod common;

use common::{benchpro_bin, unique_temp_dir};

/// Directory holding the committed fixture.
fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/test_data/align")
}

fn fixture(name: &str) -> String {
    fixture_dir().join(name).to_string_lossy().into_owned()
}

/// Writes a samplesheet and runs `benchpro align` over it, returning the output prefix.
fn run_align(test_name: &str, meta_rows: &str, extra_args: &[&str]) -> PathBuf {
    let dir = unique_temp_dir(&format!("benchpro_align_{test_name}"));
    fs::create_dir_all(&dir).expect("create temp dir");

    let meta = dir.join("meta.tsv");
    fs::write(&meta, meta_rows).expect("write meta");

    let prefix = dir.join("out");
    let mut args = vec![
        "--log-level".to_string(),
        "error".to_string(),
        "align".to_string(),
        "--meta".to_string(),
        meta.to_string_lossy().into_owned(),
        "--outprefix".to_string(),
        prefix.to_string_lossy().into_owned(),
    ];
    args.extend(extra_args.iter().map(|a| a.to_string()));

    let output = Command::new(benchpro_bin())
        .args(&args)
        .output()
        .expect("run benchpro align");
    assert!(
        output.status.success(),
        "benchpro align failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    prefix
}

/// Reads a TSV into rows keyed by column name.
fn read_tsv(path: &PathBuf) -> Vec<HashMap<String, String>> {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

/// The summary row of one tool.
fn row_for<'a>(rows: &'a [HashMap<String, String>], tool: &str) -> &'a HashMap<String, String> {
    rows.iter()
        .find(|r| r["tool"] == tool)
        .unwrap_or_else(|| panic!("no summary row for {tool}"))
}

fn number(row: &HashMap<String, String>, column: &str) -> f64 {
    row.get(column)
        .unwrap_or_else(|| panic!("no column {column}"))
        .parse()
        .unwrap_or_else(|e| panic!("{column} = {:?}: {e}", row.get(column)))
}

/// A samplesheet naming both SAM tools under `full` scoring.
fn full_scoring_meta() -> String {
    format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tReference\tScoring\n\
         ds\ts1\tgood\t{good}\t{truth}\t{c2g}\t{reference}\tfull\n\
         ds\ts1\tsloppy\t{sloppy}\t{truth}\t{c2g}\t{reference}\tfull\n",
        good = fixture("good.sam"),
        sloppy = fixture("sloppy.sam"),
        truth = fixture("reads.truth.tsv"),
        c2g = fixture("reference.contig2genome.tsv"),
        reference = fixture("reference.fna"),
    )
}

#[test]
fn full_scoring_strata_match_the_fixture_by_construction() {
    let prefix = run_align("full", &full_scoring_meta(), &["--verify-sample", "0"]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));

    let good = row_for(&summary, "good");
    // 8 pairs placed exactly + 1 pair on the sibling marker = 18 of 20 mates; pair 10 is absent.
    assert_eq!(number(good, "total"), 20.0);
    assert_eq!(number(good, "aligned"), 18.0);
    // Every placement is on genomeA, including the one on the wrong marker of that genome.
    assert_eq!(number(good, "correct"), 18.0);
    assert_eq!(number(good, "align_pct"), 90.0);
    assert_eq!(number(good, "correct_pct"), 100.0);
    // Recall divides by ALL reads in the truth, not by what the tool placed.
    assert_eq!(number(good, "recall_pct"), 90.0);
    // Only 16 of the 20 truth mates are on the true contig: pair 9 went to 1001_geneB and pair 10
    // was not placed. The strata are shares of the truth, not of what the tool placed.
    assert!((number(good, "reference_pct") - 100.0 * 16.0 / 20.0).abs() < 1e-9);
    assert!((number(good, "position_pct") - 100.0 * 16.0 / 20.0).abs() < 1e-9);
    // Its precision counterpart divides by the 18 it did place.
    assert!((number(good, "position_precision_pct") - 100.0 * 16.0 / 18.0).abs() < 1e-9);

    let sloppy = row_for(&summary, "sloppy");
    // 10 (pairs 1-5) + 6 (pairs 6-8, wrong genome) + 2 (pair 9) + 1 (pair 10 mate 2) = 19.
    // The unmapped mate and the secondary record are not placements.
    assert_eq!(number(sloppy, "aligned"), 19.0);
    assert_eq!(number(sloppy, "correct"), 13.0);
    assert_eq!(number(sloppy, "align_pct"), 95.0);
    assert!((number(sloppy, "correct_pct") - 100.0 * 13.0 / 19.0).abs() < 1e-9);
    assert_eq!(number(sloppy, "aln_unmapped"), 1.0);
    assert_eq!(number(sloppy, "aln_secondary"), 1.0);
}

#[test]
fn species_scoring_ignores_which_marker_of_the_species_was_hit() {
    // Same fixture, but truth's genome column is compared against the contig prefix. Under species
    // scoring the truth genome must be the species id, so a dedicated truth is not needed here:
    // the fixture's genome label is `genomeA` and the contig prefix is `1001`, which never match --
    // so this asserts the *opposite* direction: species scoring on a contig2genome fixture is all
    // wrong, and the module says so rather than silently agreeing with the full-scoring answer.
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tScoring\n\
         ds\ts1\tgood\t{good}\t{truth}\tspecies\n",
        good = fixture("good.sam"),
        truth = fixture("reads.truth.tsv"),
    );
    let prefix = run_align("species", &meta, &[]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let good = row_for(&summary, "good");

    assert_eq!(
        number(good, "aligned"),
        18.0,
        "placement is mode independent"
    );
    assert_eq!(
        number(good, "correct"),
        0.0,
        "prefix '1001' is not the label 'genomeA'"
    );
    // The stratum columns exist only under full scoring.
    assert!(!good.contains_key("position_pct"));
}

#[test]
fn paf_contenders_are_scored_at_the_mapping_level() {
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tScoring\n\
         ds\ts1\tmapper\t{paf}\t{truth}\t{c2g}\tfull\n",
        paf = fixture("mapper.paf"),
        truth = fixture("reads.truth.tsv"),
        c2g = fixture("reference.contig2genome.tsv"),
    );
    let prefix = run_align("paf", &meta, &[]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let mapper = row_for(&summary, "mapper");

    // Pairs 1-9 at the true locus; the /1,/2 suffix must join with the truth's mate column.
    assert_eq!(number(mapper, "aligned"), 18.0);
    assert_eq!(number(mapper, "correct"), 18.0);
    // 18 of the 20 truth mates are at their true locus -- pair 10 was not placed at all.
    assert_eq!(number(mapper, "position_pct"), 90.0);
    // Of the 18 it did place, every one is right.
    assert_eq!(number(mapper, "position_precision_pct"), 100.0);
    // PAF carries no CIGAR, SEQ or NM, so there is nothing base level to report.
    assert_eq!(mapper["aln_identity"], "");
    assert_eq!(mapper["aln_verified"], "");
}

#[test]
fn replay_recomputes_identity_and_catches_a_lying_nm_tag() {
    let prefix = run_align("replay", &full_scoring_meta(), &["--verify-sample", "0"]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));

    let good = row_for(&summary, "good");
    // Every retained record is replayed.
    assert_eq!(number(good, "aln_verified"), 18.0);
    // Pairs 1-8 are exact copies of the reference, so they replay to identity 100%. Pair 9 sits on
    // the wrong marker, whose bases differ, so the mean is below 100%.
    assert!(number(good, "aln_identity_v") > 80.0);
    assert!(number(good, "aln_identity_v") < 100.0);
    // The tool claims NM:i:0 on pairs 1-8 (true) and NM:i:3 on pair 9 (almost certainly not the
    // real distance to the other marker), so agreement is high but not total.
    assert!(number(good, "aln_nm_agree") >= 100.0 * 16.0 / 18.0 - 1e-9);
    assert!(number(good, "aln_nm_agree") < 100.0);

    let sloppy = row_for(&summary, "sloppy");
    // Pairs 1-5 emit no NM at all: 10 of 19 placements.
    assert!((number(sloppy, "aln_no_nm") - 100.0 * 10.0 / 19.0).abs() < 1e-9);
    // One malformed record of 19.
    assert!((number(sloppy, "aln_malformed") - 100.0 / 19.0).abs() < 1e-9);
    // Reported identity covers only the 9 records that carry an NM tag; a replayed identity covers
    // all of them, which is exactly the point of the replay.
    assert!(number(sloppy, "aln_verified") > 0.0);
}

#[test]
fn the_head_to_head_compares_the_two_contenders() {
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tReference\tScoring\tPeer\n\
         ds\ts1\tgood\t{good}\t{truth}\t{c2g}\t{reference}\tfull\tsloppy\n\
         ds\ts1\tsloppy\t{sloppy}\t{truth}\t{c2g}\t{reference}\tfull\t\n",
        good = fixture("good.sam"),
        sloppy = fixture("sloppy.sam"),
        truth = fixture("reads.truth.tsv"),
        c2g = fixture("reference.contig2genome.tsv"),
        reference = fixture("reference.fna"),
    );
    let prefix = run_align("h2h", &meta, &["--verify-sample", "0"]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));

    let good = row_for(&summary, "good");
    assert_eq!(good["h2h_peer"], "sloppy");
    // good placed 18 mates, sloppy 19; they overlap on pairs 1-9 minus what each missed.
    assert!(number(good, "h2h_common") > 0.0);
    assert!(
        number(good, "h2h_only_peer") > 0.0,
        "sloppy placed pair 10 mate 2"
    );
    // Only the peer's row lacks a head-to-head, because only `good` named a Peer.
    let sloppy = row_for(&summary, "sloppy");
    assert_eq!(sloppy["h2h_peer"], "");
}

#[test]
fn the_mapq_curve_trades_recall_for_precision() {
    let prefix = run_align("mapq", &full_scoring_meta(), &[]);
    let curve = read_tsv(&prefix.with_extension("align_mapq.tsv"));

    let good: Vec<&HashMap<String, String>> =
        curve.iter().filter(|r| r["tool"] == "good").collect();
    assert!(!good.is_empty());

    // Cutoffs ascend, and `kept` can only shrink as the cutoff rises.
    let mut previous_kept = f64::INFINITY;
    for point in &good {
        let kept = number(point, "kept");
        assert!(
            kept <= previous_kept,
            "kept must be monotonically non-increasing"
        );
        previous_kept = kept;
    }
    // good places pairs 1-8 at MAPQ 60 and pair 9 at 40, so filtering at 60 drops the wrong-marker
    // pair and precision reaches 100%.
    let top = good.last().expect("a curve point");
    assert_eq!(number(top, "mapq"), 60.0);
    assert_eq!(number(top, "precision_pct"), 100.0);

    // Both recall denominators are reported, and they differ: the field-relative one divides by
    // what the best contender achieved, the absolute one by the whole truth.
    assert!(number(top, "recall_mappable_pct") >= number(top, "recall_total_pct"));
    for column in [
        "recall_mappable_pct",
        "recall_total_pct",
        "f1_mappable",
        "f1_total",
    ] {
        let value = number(top, column);
        assert!((0.0..=100.0).contains(&value), "{column} = {value}");
    }
}

#[test]
fn per_read_output_covers_every_truth_read_for_every_tool() {
    let prefix = run_align("per_read", &full_scoring_meta(), &["--per-read"]);
    let reads = read_tsv(&prefix.with_extension("align_reads.tsv"));

    assert_eq!(reads.len(), 40, "20 truth mates x 2 tools");

    let good: Vec<&HashMap<String, String>> =
        reads.iter().filter(|r| r["tool"] == "good").collect();
    let count = |verdict: &str| good.iter().filter(|r| r["verdict"] == verdict).count();

    assert_eq!(count("position"), 16);
    assert_eq!(count("genome"), 2, "pair 9, on the sibling marker");
    assert_eq!(count("unaligned"), 2, "pair 10");
    assert_eq!(count("reference"), 0);
}

#[test]
fn base_level_columns_are_shares_of_emitted_alignments_not_of_truth_reads() {
    // A tool routinely places reads the truth does not cover: a subsampled truth, a marker DB, a
    // spike-in. `aln_*` must then divide by what the tool emitted (19 primaries here), not by the
    // 4 truth reads it happened to place -- dividing by the latter yields percentages over 100.
    let dir = unique_temp_dir("benchpro_align_denominator");
    fs::create_dir_all(&dir).expect("create temp dir");
    let truth = dir.join("four.truth.tsv");
    let full = fs::read_to_string(fixture("reads.truth.tsv")).expect("read truth");
    let head: String = full.lines().take(4).map(|l| format!("{l}\n")).collect();
    fs::write(&truth, head).expect("write truth");

    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tScoring\n\
         ds\ts1\tsloppy\t{sloppy}\t{truth}\tspecies\n",
        sloppy = fixture("sloppy.sam"),
        truth = truth.to_string_lossy(),
    );
    let prefix = run_align("denominator", &meta, &[]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let row = row_for(&summary, "sloppy");

    assert_eq!(
        number(row, "aligned"),
        4.0,
        "only 4 truth reads were placed"
    );
    assert_eq!(
        number(row, "aln_records"),
        19.0,
        "but 19 alignments were emitted"
    );
    // 10 of the 19 primaries carry no NM tag.
    assert!((number(row, "aln_no_nm") - 100.0 * 10.0 / 19.0).abs() < 1e-9);

    for column in [
        "aln_no_nm",
        "aln_proper_pair",
        "aln_malformed",
        "aln_coverage",
    ] {
        let value = number(row, column);
        assert!(
            (0.0..=100.0).contains(&value),
            "{column} = {value} is not a percentage"
        );
    }
}

#[test]
fn the_per_read_table_covers_exactly_the_samples_the_summary_kept() {
    // `good` ran on two samples, `sloppy` on one. s2 is dropped from the aggregates because not
    // every contender produced it -- so the per-read table must drop it too, or it describes runs
    // that no summary row accounts for.
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tScoring\n\
         ds\ts1\tgood\t{good}\t{truth}\tspecies\n\
         ds\ts2\tgood\t{good}\t{truth}\tspecies\n\
         ds\ts1\tsloppy\t{sloppy}\t{truth}\tspecies\n",
        good = fixture("good.sam"),
        sloppy = fixture("sloppy.sam"),
        truth = fixture("reads.truth.tsv"),
    );
    let prefix = run_align("per_read_restricted", &meta, &["--per-read"]);

    let samples = read_tsv(&prefix.with_extension("align_samples.tsv"));
    assert!(
        samples.iter().all(|r| r["sample"] == "s1"),
        "s2 should have been dropped from the aggregates"
    );

    let reads = read_tsv(&prefix.with_extension("align_reads.tsv"));
    assert!(!reads.is_empty());
    assert!(
        reads.iter().all(|r| r["sample"] == "s1"),
        "per-read rows survive for a sample the summary dropped"
    );
}

#[test]
fn a_gold_standard_sam_answers_all_three_benchmarks() {
    // Read simulators (art_illumina -sam) emit a SAM recording not just where each read came from
    // but how it truly aligns. That third fact is what a truth TSV cannot express, so it unlocks
    // the third question: not "near the right place" but "the alignment the simulator produced".
    //
    // approx.sam, against gold.sam, is built so the three answers differ:
    //   pairs 1-6  identical to gold                    -> exact       (12 mates)
    //   pairs 7-8  right locus, different CIGAR         -> position    (4 mates)
    //   pair  9    right contig, 30 bp away             -> reference   (2 mates)
    //   pair  10   wrong genome                         -> wrong       (2 mates)
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tScoring\n\
         ds\ts1\tapprox\t{approx}\t{gold}\t{c2g}\tfull\n",
        approx = fixture("approx.sam"),
        gold = fixture("gold.sam"),
        c2g = fixture("reference.contig2genome.tsv"),
    );
    let prefix = run_align("gold_sam", &meta, &[]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let row = row_for(&summary, "approx");

    assert_eq!(number(row, "total"), 20.0);
    assert_eq!(number(row, "aligned"), 20.0);

    // (a) the right reference: everything except pair 10, which went to the other genome.
    assert_eq!(number(row, "reference_pct"), 100.0 * 18.0 / 20.0);
    // (b) the right position within the window: pair 9 is 30 bp out, beyond the default 100? no --
    //     30 <= 100, so it counts. Pairs 1-8 and 9 = 18 mates.
    assert_eq!(number(row, "position_pct"), 100.0 * 18.0 / 20.0);
    // (c) the identical alignment: only the 12 mates of pairs 1-6.
    assert_eq!(number(row, "exact_pct"), 100.0 * 12.0 / 20.0);
}

#[test]
fn the_exact_benchmark_separates_a_shifted_alignment_from_a_reshaped_one() {
    // Tightening the window to 0 makes (b) demand the exact start, which pair 9 (30 bp out) fails.
    // (c) is still stricter: pairs 7-8 share gold's start but not its CIGAR.
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tScoring\n\
         ds\ts1\tapprox\t{approx}\t{gold}\t{c2g}\tfull\n",
        approx = fixture("approx.sam"),
        gold = fixture("gold.sam"),
        c2g = fixture("reference.contig2genome.tsv"),
    );
    let prefix = run_align("gold_sam_tight", &meta, &["--tolerance", "0"]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let row = row_for(&summary, "approx");

    // pair 9 now fails the window, leaving pairs 1-8 = 16 mates.
    assert_eq!(number(row, "position_pct"), 100.0 * 16.0 / 20.0);
    // ...and of those, only pairs 1-6 reproduce gold's CIGAR.
    assert_eq!(number(row, "exact_pct"), 100.0 * 12.0 / 20.0);
    assert!(
        number(row, "exact_pct") < number(row, "position_pct"),
        "the exact benchmark must be strictly harder than the windowed one"
    );
}

#[test]
fn a_truth_tsv_leaves_the_exact_benchmark_undefined() {
    // A TSV records where a read came from, not how it aligns, so the third column must stay empty
    // rather than reporting 0% -- which would read as "never reproduces the gold alignment".
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tScoring\n\
         ds\ts1\tgood\t{good}\t{truth}\t{c2g}\tfull\n",
        good = fixture("good.sam"),
        truth = fixture("reads.truth.tsv"),
        c2g = fixture("reference.contig2genome.tsv"),
    );
    let prefix = run_align("tsv_no_exact", &meta, &[]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let row = row_for(&summary, "good");

    assert_eq!(
        row["exact_pct"], "",
        "a TSV truth cannot answer benchmark (c)"
    );
    assert!(
        number(row, "position_pct") > 0.0,
        "but it still answers (a) and (b)"
    );
}

#[test]
fn a_contig_missing_from_the_map_is_not_scored_as_the_wrong_genome() {
    // The scorer labels an unmapped contig "NA"; the gold-SAM loader must use the same word, or a
    // reference with plasmids or decoys -- an incomplete map is the normal state -- scores 0%.
    let dir = unique_temp_dir("benchpro_align_partial_map");
    fs::create_dir_all(&dir).expect("create temp dir");
    let map = dir.join("partial.tsv");
    fs::write(&map, "2002_geneA\tgenomeB\n").expect("write map");

    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tScoring\n\
         ds\ts1\tself\t{gold}\t{gold}\t{map}\tfull\n",
        gold = fixture("gold.sam"),
        map = map.to_string_lossy(),
    );
    let prefix = run_align("partial_map", &meta, &[]);
    let row_data = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let row = row_for(&row_data, "self");

    assert_eq!(
        number(row, "correct_pct"),
        100.0,
        "a file must score 100% against itself"
    );
    assert_eq!(number(row, "reference_pct"), 100.0);
    assert_eq!(number(row, "exact_pct"), 100.0);
}

#[test]
fn species_scoring_labels_a_gold_sam_with_the_species_prefix() {
    // `species` never consults a contig map, so the truth label has to be written in the vocabulary
    // the scorer uses: the contig's prefix. Labelling it with the full contig name scores 0%.
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tScoring\n\
         ds\ts1\tself\t{gold}\t{gold}\tspecies\n",
        gold = fixture("gold.sam"),
    );
    let prefix = run_align("species_gold", &meta, &[]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let row = row_for(&summary, "self");

    assert_eq!(number(row, "correct_pct"), 100.0);
    assert_eq!(number(row, "recall_pct"), 100.0);
}

#[test]
fn a_paf_truth_leaves_the_exact_benchmark_undefined() {
    // PAF carries no CIGAR, so "identical to the gold alignment" is unanswerable. Reporting 100%
    // (neither side has one) or 0% (the tool never matches) would both be inventions.
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tScoring\n\
         ds\ts1\tp\t{paf}\t{paf}\t{c2g}\tfull\n",
        paf = fixture("mapper.paf"),
        c2g = fixture("reference.contig2genome.tsv"),
    );
    let prefix = run_align("paf_truth", &meta, &[]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let row = row_for(&summary, "p");

    assert_eq!(
        row["exact_pct"], "",
        "a PAF truth cannot answer benchmark (c)"
    );
    assert_eq!(
        number(row, "position_pct"),
        100.0,
        "but it answers (a) and (b)"
    );
}

#[test]
fn the_per_read_table_has_a_header_even_when_no_sample_survives() {
    // Two tools with disjoint samples: the intersection is empty, so nothing is scored. Every
    // other table still gets its header, and a zero-byte file is not a table -- readers that
    // return an empty frame for the others raise on it.
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tScoring\n\
         ds\ts1\ta\t{good}\t{truth}\tspecies\n\
         ds\ts2\tb\t{good}\t{truth}\tspecies\n",
        good = fixture("good.sam"),
        truth = fixture("reads.truth.tsv"),
    );
    let prefix = run_align("empty_per_read", &meta, &["--per-read"]);
    let reads = prefix.with_extension("align_reads.tsv");

    let text = fs::read_to_string(&reads).expect("read per-read table");
    assert!(!text.is_empty(), "an empty file is not a table");
    assert!(
        text.starts_with("dataset\tsample\ttool\t"),
        "header missing: {text:?}"
    );
    assert_eq!(text.lines().count(), 1, "header only, no rows");
}

#[test]
fn indel_recovery_is_reported_separately_from_the_easy_majority() {
    // approx.sam reproduces pairs 1-6 exactly and reshapes pairs 7-8. In gold.sam pair 9 carries a
    // deletion and pair 10 a soft clip, so the hard subsets are small and distinct from the total.
    let meta = format!(
        "ID\tSample\tTool\tAlignment\tTruth\tContig2Genome\tScoring\n\
         ds\ts1\tapprox\t{approx}\t{gold}\t{c2g}\tfull\n",
        approx = fixture("approx.sam"),
        gold = fixture("gold.sam"),
        c2g = fixture("reference.contig2genome.tsv"),
    );
    let prefix = run_align("indel_recovery", &meta, &[]);
    let summary = read_tsv(&prefix.with_extension("align_summary.tsv"));
    let row = row_for(&summary, "approx");

    // Pair 9 = 2 mates carry a deletion in the gold alignment; pair 10 = 2 mates are clipped.
    assert_eq!(number(row, "indel_reads"), 2.0);
    assert_eq!(number(row, "clipped_reads"), 2.0);
    // approx puts pair 9 thirty bases out with a plain 20M, so it never reproduces the deletion.
    assert_eq!(number(row, "indel_exact_pct"), 0.0);
    // ...and pair 10 it places on the wrong genome entirely.
    assert_eq!(number(row, "clipped_position_pct"), 0.0);
    // The overall exact rate is far higher, which is exactly why the subsets are reported.
    assert!(number(row, "exact_pct") > 50.0);
}

#[test]
fn a_truth_that_shares_no_read_is_an_error_not_a_zero_score() {
    let dir = unique_temp_dir("benchpro_align_mismatch");
    fs::create_dir_all(&dir).expect("create temp dir");
    let truth = dir.join("other.truth.tsv");
    fs::write(&truth, "totally_different_id\t1\tctg\t1\tgenomeA\n").expect("write truth");

    let meta = dir.join("meta.tsv");
    fs::write(
        &meta,
        format!(
            "ID\tSample\tTool\tAlignment\tTruth\tScoring\n\
             ds\ts1\tgood\t{}\t{}\tspecies\n",
            fixture("good.sam"),
            truth.to_string_lossy()
        ),
    )
    .expect("write meta");

    let output = Command::new(benchpro_bin())
        .args([
            "--log-level",
            "error",
            "align",
            "--meta",
            &meta.to_string_lossy(),
            "--outprefix",
            &dir.join("out").to_string_lossy(),
        ])
        .output()
        .expect("run benchpro align");

    assert!(
        !output.status.success(),
        "a zero overlap must not look like a clean run"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("share no read id"), "stderr was: {stderr}");
}

#[test]
fn an_invalid_samplesheet_reports_every_problem_at_once() {
    let dir = unique_temp_dir("benchpro_align_badmeta");
    fs::create_dir_all(&dir).expect("create temp dir");
    let meta = dir.join("meta.tsv");
    fs::write(
        &meta,
        "ID\tSample\tTool\tAlignment\tTruth\tPeer\n\
         ds\ts1\tgood\t/no/such/file.sam\t/no/such/truth.tsv\tghost\n",
    )
    .expect("write meta");

    let output = Command::new(benchpro_bin())
        .args([
            "--log-level",
            "error",
            "align",
            "--meta",
            &meta.to_string_lossy(),
            "--outprefix",
            &dir.join("out").to_string_lossy(),
        ])
        .output()
        .expect("run benchpro align");

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    // All three problems, not just the first one found.
    assert!(stderr.contains("file.sam' does not exist"), "{stderr}");
    assert!(stderr.contains("truth.tsv' does not exist"), "{stderr}");
    assert!(stderr.contains("Peer 'ghost'"), "{stderr}");
}

#[test]
fn validate_meta_checks_without_scoring_anything() {
    let prefix = run_align("validate", &full_scoring_meta(), &["--validate-meta"]);
    assert!(
        !prefix.with_extension("align_summary.tsv").exists(),
        "--validate-meta must not write output tables"
    );
}
