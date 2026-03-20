use std::{
    fs,
    path::PathBuf,
    process::Command,
};

mod common;

use common::{benchpro_bin, unique_temp_dir};

#[test]
fn test_profile_output_matches_golden() {
    let bin = benchpro_bin();

    let temp_dir = unique_temp_dir("benchpro_profile_test");
    fs::create_dir_all(&temp_dir).expect("Failed to create temp output dir");

    let outprefix = temp_dir.join("test");

    let status = Command::new(&bin)
        .args([
            "--log-level",
            "error",
            "profile",
            "--meta",
            "data/test_data/meta/meta.xlsx",
            "--outprefix",
            outprefix.to_str().expect("Non-UTF8 outprefix"),
            "--adjusted",
            "--allow-alternatives",
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
        summary_text.contains("MG-TK-v2_gastrointestinal_0"),
        "Expected MG-TK sample not found in summary output"
    );
    assert!(
        summary_text.contains("sylph_gastrointestinal_4"),
        "Expected sylph sample not found in summary output"
    );
    assert!(
        detailed_text.contains("MG-TK-v2_gastrointestinal_0"),
        "Expected MG-TK sample not found in detailed output"
    );
}
