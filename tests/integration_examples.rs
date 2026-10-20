//! Runs every `example-*` recipe of the Justfile and checks the lines its README promises.
//! Needs `just` on PATH (installed in CI).

use std::{fs, path::Path, process::Command};

mod common;

use common::{benchpro_bin, unique_temp_dir};

/// Runs the `just` recipe for `examples/<demo>` and returns its combined stdout and stderr.
fn run_demo(demo: &str) -> String {
    let recipe = format!("example-{}", demo.split_once('-').unwrap().1);
    // An empty HOME keeps a local ~/.taxdet/index from changing the output.
    let home = unique_temp_dir(&format!("benchpro_examples_home_{}", demo));
    fs::create_dir_all(&home).unwrap();
    let output = Command::new("just")
        .arg(&recipe)
        .env("BENCHPRO", benchpro_bin())
        .env("HOME", &home)
        .output()
        .expect("Failed to run `just`; install it with `cargo install just`");
    fs::remove_dir_all(&home).unwrap();
    let text = format!(
        "{}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.status.success(), "{} failed:\n{}", demo, text);
    text
}

fn assert_contains(demo: &str, text: &str, expected: &[&str]) {
    for line in expected {
        assert!(
            text.contains(line),
            "{}: missing '{}' in:\n{}",
            demo,
            line,
            text
        );
    }
}

fn assert_files(demo: &str, files: &[&str]) {
    for file in files {
        let path = Path::new("examples").join(demo).join("output").join(file);
        assert!(path.is_file(), "{}: missing {}", demo, path.display());
    }
}

#[test]
fn test_example_01_minimal() {
    let text = run_demo("01-minimal");
    assert_contains(
        "01-minimal",
        &text,
        &[
            "Row 'sylph_3': Taxonomy GTDB (inferred)",
            "Done: 4 rows, taxonomy given for 0, inferred for 4",
            "Warnings: 0",
        ],
    );
    assert_files(
        "01-minimal",
        &["run.tsv", "run_detailed.tsv", "run_inputs.tsv"],
    );
}

#[test]
fn test_example_02_full() {
    let text = run_demo("02-full");
    assert_contains(
        "02-full",
        &text,
        &[
            "Adjusted: 2 rows, skipped 0 (no GoldStdTree)",
            "Warnings: 0",
        ],
    );
    let detailed = fs::read_to_string("examples/02-full/output/run_detailed.tsv").unwrap();
    assert_eq!(
        detailed.matches("\tFFP\t").count(),
        3,
        "expected 3 FFP rows"
    );
}

#[test]
fn test_example_03_gtdb_and_ncbi() {
    let text = run_demo("03-gtdb-and-ncbi");
    assert_contains(
        "03-gtdb-and-ncbi",
        &text,
        &[
            "Taxonomy inferred as GTDB for 2 row(s): MetaPhlAn4_GTDB_4, sylph_4",
            "Taxonomy inferred as NCBI for 2 row(s): MetaPhlAn4_4, mOTUs3_4",
            "Warnings: 0",
        ],
    );
}

#[test]
fn test_example_04_mistakes() {
    let text = run_demo("04-mistakes");
    assert_contains(
        "04-mistakes",
        &text,
        &[
            "missing required column(s): GoldStd",
            "prediction looks like GTDB but gold standard looks like NCBI",
            "tool 'MetaPhlAn4_GTDB' mixes taxonomies",
            "row 2: column 'Tool' is empty",
            "--adjusted needs a GoldStdTree for at least one row",
            "Output folder 'examples/04-mistakes/no-such-folder' does not exist",
        ],
    );
    assert_eq!(text.matches("--> exit code 1").count(), 4, "{}", text);
    assert_eq!(text.matches("--> exit code 2").count(), 2, "{}", text);
    assert!(!text.contains("panicked"), "{}", text);
}

#[test]
fn test_example_05_matrix() {
    let text = run_demo("05-matrix");
    assert_contains(
        "05-matrix",
        &text,
        &["Done: 20 rows, taxonomy given for 20, inferred for 0"],
    );
}

#[test]
fn test_example_06_other_commands() {
    let text = run_demo("06-other-commands");
    assert_contains(
        "06-other-commands",
        &text,
        &["tool\tsylph\tX", "tool\tmetaphlan\t4", "tool\tmotus\t4.0.4"],
    );
    assert_files(
        "06-other-commands",
        &[
            "sylph.cami.profile",
            "sylph_species_matrix.tsv",
            "strain.monophyly.tsv",
            "strain.msa_error.tsv",
        ],
    );
}
