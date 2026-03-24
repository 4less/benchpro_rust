use std::{fs, process::Command};

mod common;

use common::{benchpro_bin, unique_temp_dir};

fn run_detect(bin: &std::path::Path, input: &str) -> String {
    let output = Command::new(bin)
        .args(["normalize", "--input", input, "--detect"])
        .output()
        .expect("Failed to run benchpro normalize --detect");
    assert!(
        output.status.success(),
        "normalize --detect failed for {input}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .expect("stdout was not valid UTF-8")
        .trim()
        .to_owned()
}

fn run_normalize_plain(bin: &std::path::Path, input: &str, output: &std::path::Path) -> String {
    let status = Command::new(bin)
        .args([
            "normalize",
            "--input",
            input,
            "--output",
            output.to_str().expect("Non-UTF8 output path"),
        ])
        .status()
        .expect("Failed to run benchpro normalize");
    assert!(status.success(), "normalize failed for {input}");
    fs::read_to_string(output).expect("Failed to read normalized output")
}

#[test]
fn test_normalize_detect_regression_examples() {
    let bin = benchpro_bin();
    let cases = [
        (
            "data/profile_examples/metaphlan3013_camioutput.profile",
            "cami\tmetaphlan\t0.10.0",
        ),
        (
            "data/profile_examples/metaphlan402.profile",
            "tool\tmetaphlan\t4",
        ),
        ("data/profile_examples/sylph090.profile", "tool\tsylph\tX"),
        ("data/profile_examples/protal.profile", "tool\tprotal\tX"),
        ("data/profile_examples/motus404.relab", "tool\tmotus\t4.0.4"),
    ];

    for (input, expected) in cases {
        let observed = run_detect(&bin, input);
        assert_eq!(observed, expected, "detection mismatch for {input}");
    }
}

#[test]
fn test_normalize_plain_regression_examples() {
    let bin = benchpro_bin();
    let temp_dir = unique_temp_dir("benchpro_normalize_regression");
    fs::create_dir_all(&temp_dir).expect("Failed to create temp output dir");

    let cases = [
        "data/profile_examples/metaphlan402.profile",
        "data/profile_examples/sylph090.profile",
        "data/profile_examples/protal.profile",
        "data/profile_examples/motus404.relab",
        "data/profile_examples/metaphlan3013_camioutput.profile",
    ];

    for input in cases {
        let safe_name = input
            .rsplit('/')
            .next()
            .unwrap_or("out")
            .replace('.', "_")
            .replace('-', "_");
        let out = temp_dir.join(format!("{safe_name}.norm.tsv"));
        let text = run_normalize_plain(&bin, input, &out);
        let mut lines = text.lines();

        let source = lines.next().unwrap_or_default();
        let format = lines.next().unwrap_or_default();
        let tool = lines.next().unwrap_or_default();
        let version = lines.next().unwrap_or_default();
        let header = lines.next().unwrap_or_default();

        assert!(
            source.starts_with("#source_profile\t"),
            "missing source header"
        );
        assert!(
            source.ends_with(input),
            "source path header should match input for {input}"
        );
        assert!(
            format.starts_with("#detected_format\t"),
            "missing format header"
        );
        assert!(tool.starts_with("#detected_tool\t"), "missing tool header");
        assert!(
            version.starts_with("#detected_version\t"),
            "missing version header"
        );
        assert_eq!(
            header, "identifier\tlineage\tabundance\tvertical_coverage",
            "unexpected data header for {input}"
        );

        let mut saw_entry = false;
        for row in lines {
            if row.trim().is_empty() {
                continue;
            }
            saw_entry = true;
            let cols = row.split('\t').collect::<Vec<_>>();
            assert!(cols.len() >= 4, "expected >=4 columns for {input}");
            assert!(
                !cols[0].trim().is_empty(),
                "identifier is empty in normalized row for {input}"
            );
            assert!(
                !cols[1].trim().is_empty(),
                "lineage is empty in normalized row for {input}"
            );
            assert!(
                cols[2].trim().parse::<f64>().is_ok(),
                "abundance is not numeric for {input}: {}",
                cols[2]
            );
            let lineage_tokens = cols[1]
                .split(['|', ';'])
                .map(str::trim)
                .filter(|token| !token.is_empty())
                .collect::<Vec<_>>();
            assert!(
                !lineage_tokens.is_empty(),
                "lineage has no tokens in normalized row for {input}"
            );
            for token in lineage_tokens {
                if let Some((prefix, _)) = token.split_once("__") {
                    assert!(
                        matches!(prefix, "k" | "d" | "p" | "c" | "o" | "f" | "g" | "s" | "t"),
                        "non-standard lineage rank prefix '{prefix}' for {input}"
                    );
                }
            }
        }
        assert!(saw_entry, "no normalized entries were written for {input}");
    }
}
