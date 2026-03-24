use std::{fs, io::Write, path::PathBuf, process::Command};

mod common;

use common::{benchpro_bin, unique_temp_dir};

#[test]
fn test_matrix_output_matches_golden() {
    let bin = benchpro_bin();
    let repo_root = std::env::current_dir().expect("Failed to get current directory");

    let temp_dir = unique_temp_dir("benchpro_matrix_test");
    fs::create_dir_all(&temp_dir).expect("Failed to create temp output dir");
    let meta_path = temp_dir.join("meta_matrix.csv");

    let outprefix = temp_dir.join("test_matrix");

    let profile_matrix =
        repo_root.join("data/test_data/profiles/predictions/protal/abundance_matrix.tsv");
    let gold_matrix =
        repo_root.join("data/test_data/profiles/gold_standard/GTDB/r214/abundance_matrix.tsv");
    let tree_path = repo_root.join("data/test_data/trees/gtdb/r214/bac120_r214.sp_labels.tree");
    let detectable_path =
        repo_root.join("data/test_data/detectable_species/gtdb_all_species_r214.tsv");

    let mut meta_file = fs::File::create(&meta_path).expect("Failed to create matrix meta fixture");
    writeln!(
        meta_file,
        "ID,Sample,Dataset,Tool,Taxonomy,Profile,ProfileRegex,GoldStd,GoldStdRegex,GoldStdTree,AvailableSpecies"
    )
    .expect("Failed to write matrix meta header");
    writeln!(
        meta_file,
        "protal_Gastrointestinal,NA,Gastrointestinal,protal,GTDB r214,{},^GAST([0-9]+)$,{},^gastrointestinal_([0-9]+)$,{},{}",
        profile_matrix.display(),
        gold_matrix.display(),
        tree_path.display(),
        detectable_path.display()
    )
    .expect("Failed to write matrix meta row");

    let status = Command::new(&bin)
        .args([
            "--log-level",
            "error",
            "profile",
            "--meta",
            meta_path.to_str().expect("Non-UTF8 meta path"),
            "--outprefix",
            outprefix.to_str().expect("Non-UTF8 outprefix"),
            "--adjusted",
            "--allow-alternatives",
            "--ignore-abundance-error",
        ])
        .status()
        .expect("Failed to run benchpro");

    assert!(status.success(), "benchpro exited with a non-zero status");

    let actual_summary = outprefix.with_extension("tsv");
    let actual_detailed = PathBuf::from(format!("{}_detailed.tsv", outprefix.display()));

    let summary_text = fs::read_to_string(&actual_summary).expect("Failed to read summary output");
    let detailed_text =
        fs::read_to_string(&actual_detailed).expect("Failed to read detailed output");

    assert!(
        summary_text.starts_with("Rank\tID\tAllowAlternatives\tAdjusted\t"),
        "Unexpected summary header"
    );
    assert!(
        detailed_text.starts_with("Name\tType\tPredictionAbundance\tGoldStdAbundance"),
        "Unexpected detailed header"
    );
    assert!(
        summary_text.contains("protal_Gastrointestinal_0"),
        "Expected gastrointestinal matrix sample not found in summary output"
    );
    assert!(
        detailed_text.contains("protal_Gastrointestinal_0"),
        "Expected gastrointestinal matrix sample not found in detailed output"
    );
}
