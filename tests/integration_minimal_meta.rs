use std::{collections::HashMap, fs, path::Path, process::Command};

mod common;

use common::{benchpro_bin, unique_temp_dir};

const SYLPH: &str = "data/test_data/profiles/predictions/sylph/Gastrointestinal_4.profile";
const PROTAL: &str = "data/test_data/profiles/predictions/protal/4.profile";
const GOLD: &str = "data/test_data/profiles/gold_standard/GTDB/r214/abundance_4_gtdb.profile";

fn run(meta: &Path, outprefix: &Path) {
    let status = Command::new(benchpro_bin())
        .args([
            "--log-level",
            "error",
            "profile",
            "--no-taxonomy-index",
            "--meta",
        ])
        .arg(meta)
        .arg("--outprefix")
        .arg(outprefix)
        .status()
        .expect("Failed to run benchpro");
    assert!(status.success(), "benchpro failed for {}", meta.display());
}

/// (Rank, ID) -> (TP, FP, FN, F1) from a summary TSV.
fn metrics(path: &Path) -> HashMap<(String, String), Vec<String>> {
    let text = fs::read_to_string(path).expect("Failed to read summary");
    let mut lines = text.lines();
    let header = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let idx = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let (rank, id) = (idx("Rank"), idx("ID"));
    let wanted = ["TP", "FP", "FN", "F1"].map(idx);
    lines
        .map(|line| {
            let f = line.split('\t').collect::<Vec<_>>();
            let key = (f[rank].to_string(), f[id].to_string());
            (key, wanted.iter().map(|&i| f[i].to_string()).collect())
        })
        .collect()
}

#[test]
fn test_minimal_meta_matches_full_meta() {
    let dir = unique_temp_dir("benchpro_minimal_meta");
    fs::create_dir_all(&dir).unwrap();

    let full = dir.join("meta_full.tsv");
    fs::write(
        &full,
        format!(
            "ID\tSample\tTool\tTaxonomy\tProfile\tGoldStd\tGoldStdTree\tAvailableSpecies\n\
             sylph_4\ts4\tsylph\tGTDB r214\t{SYLPH}\t{GOLD}\tNA\tNA\n\
             protal_4\ts4\tprotal\tGTDB r214\t{PROTAL}\t{GOLD}\tNA\tNA\n"
        ),
    )
    .unwrap();
    let minimal = dir.join("meta_minimal.tsv");
    fs::write(
        &minimal,
        format!(
            "ID\tSample\tTool\tProfile\tGoldStd\n\
             sylph_4\ts4\tsylph\t{SYLPH}\t{GOLD}\n\
             protal_4\ts4\tprotal\t{PROTAL}\t{GOLD}\n"
        ),
    )
    .unwrap();

    run(&full, &dir.join("full"));
    run(&minimal, &dir.join("minimal"));

    let full_metrics = metrics(&dir.join("full.tsv"));
    let minimal_metrics = metrics(&dir.join("minimal.tsv"));
    assert!(!full_metrics.is_empty());
    assert_eq!(full_metrics, minimal_metrics);

    let inputs = fs::read_to_string(dir.join("minimal_inputs.tsv")).unwrap();
    assert!(inputs.contains("\tGTDB\tinferred\t"), "{}", inputs);
    fs::remove_dir_all(&dir).unwrap();
}
