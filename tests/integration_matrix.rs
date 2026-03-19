use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

mod common;

use common::{benchpro_bin, read_file, unique_temp_dir};

#[test]
fn test_matrix_output_matches_golden() {
    let bin = benchpro_bin();

    let temp_dir = unique_temp_dir("benchpro_matrix_test");
    fs::create_dir_all(&temp_dir).expect("Failed to create temp output dir");

    let outprefix = temp_dir.join("test_matrix");

    let status = Command::new(&bin)
        .args([
            "profile",
            "--meta",
            "data/test_data/meta/meta_matrix.xlsx",
            "--outprefix",
            outprefix.to_str().expect("Non-UTF8 outprefix"),
            "--adjusted",
            "--allow-alternatives",
            "--log-level",
            "error",
            "--ignore-abundance-error",
        ])
        .status()
        .expect("Failed to run benchpro");

    assert!(status.success(), "benchpro exited with a non-zero status");

    let expected_summary = Path::new("data/golden/test_meta_matrix.tsv");
    let expected_detailed = Path::new("data/golden/test_meta_matrix_detailed.tsv");

    let actual_summary = outprefix.with_extension("tsv");
    let actual_detailed = PathBuf::from(format!("{}_detailed.tsv", outprefix.display()));

    let expected_summary_bytes = read_file(expected_summary);
    let expected_detailed_bytes = read_file(expected_detailed);

    let actual_summary_bytes = read_file(&actual_summary);
    let actual_detailed_bytes = read_file(&actual_detailed);

    assert_eq!(
        expected_summary_bytes, actual_summary_bytes,
        "Summary output does not match golden file"
    );
    assert_eq!(
        expected_detailed_bytes, actual_detailed_bytes,
        "Detailed output does not match golden file"
    );
}
