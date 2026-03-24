use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use crate::normalize_detect::detect_profile;
use crate::normalize_loader::{load_normalized_for_kind, NormalizedEntry};
use crate::options::{NormalizeArgs, NormalizeOutputFormat};

/// Run the normalize command.
///
/// # Arguments
///
/// * `args` - Arguments for the normalize command
pub fn run(args: &NormalizeArgs) {
    let mut file = File::open(&args.input).expect("Failed to open input file");
    let mut content = String::new();
    file.read_to_string(&mut content)
        .expect("Failed to read input file");

    let detection = detect_profile(&content);

    if args.detect {
        println!(
            "{}\t{}\t{}",
            detection.format.as_str(),
            detection.tool,
            detection.version
        );
        return;
    }

    let output_path = args
        .output
        .as_deref()
        .expect("Output path is required unless --detect is provided");

    let entries = load_normalized_for_kind(
        &content,
        &detection.format,
        &detection.tool,
        &detection.version,
    )
    .unwrap_or_else(|err| panic!("Failed to load profile: {err}"));

    let mut output = File::create(output_path).expect("Failed to create output file");
    writeln!(output, "#source_profile\t{}", args.input).expect("Failed to write metadata header");
    writeln!(output, "#detected_format\t{}", detection.format.as_str())
        .expect("Failed to write metadata header");
    writeln!(output, "#detected_tool\t{}", detection.tool)
        .expect("Failed to write metadata header");
    writeln!(output, "#detected_version\t{}", detection.version)
        .expect("Failed to write metadata header");
    match args.output_format {
        NormalizeOutputFormat::Plain => {
            writeln!(output, "identifier\tlineage\tabundance\tvertical_coverage")
                .expect("Failed to write header");
            for entry in entries {
                let coverage = entry
                    .vertical_coverage
                    .map(|value| value.to_string())
                    .unwrap_or_default();
                writeln!(
                    output,
                    "{}\t{}\t{}\t{}",
                    entry.identifier, entry.lineage, entry.abundance, coverage
                )
                .expect("Failed to write entry");
            }
        }
        NormalizeOutputFormat::Cami => {
            write_cami(&mut output, &args.input, &detection.version, &entries);
        }
    }
}

fn write_cami<W: Write>(
    output: &mut W,
    input_path: &str,
    version: &str,
    entries: &[NormalizedEntry],
) {
    let sample_id = Path::new(input_path)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("sample");

    writeln!(output, "# Taxonomic Profiling Output").expect("Failed to write CAMI header");
    writeln!(output, "@SampleID:{sample_id}").expect("Failed to write CAMI header");
    writeln!(output, "@Version:{version}").expect("Failed to write CAMI header");
    writeln!(
        output,
        "@Ranks:superkingdom|phylum|class|order|family|genus|species|strain"
    )
    .expect("Failed to write CAMI header");
    writeln!(output, "@TaxonomyID:na").expect("Failed to write CAMI header");
    writeln!(output, "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE")
        .expect("Failed to write CAMI header");

    for entry in entries {
        let rank = infer_rank_from_lineage(&entry.lineage);
        let taxpath = entry
            .metadata
            .get("taxpath_ids")
            .map(String::as_str)
            .filter(|value| !value.is_empty())
            .unwrap_or(&entry.lineage);
        let taxid = taxpath
            .rsplit('|')
            .find(|token| !token.trim().is_empty())
            .map(str::trim)
            .unwrap_or(&entry.identifier);
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}",
            taxid, rank, taxpath, entry.lineage, entry.abundance
        )
        .expect("Failed to write CAMI entry");
    }
}

fn infer_rank_from_lineage(lineage: &str) -> &'static str {
    let last = lineage
        .rsplit(['|', ';'])
        .next()
        .map(str::trim)
        .unwrap_or_default();
    if let Some((prefix, _)) = last.split_once("__") {
        return match prefix {
            "k" | "d" | "sk" => "superkingdom",
            "p" => "phylum",
            "c" => "class",
            "o" => "order",
            "f" => "family",
            "g" => "genus",
            "s" => "species",
            "t" => "strain",
            _ => "no_rank",
        };
    }
    "no_rank"
}

#[cfg(test)]
mod tests {
    use super::{infer_rank_from_lineage, write_cami};
    use crate::normalize_loader::NormalizedEntry;
    use std::collections::BTreeMap;

    #[test]
    fn infers_strain_rank_for_t_prefix() {
        assert_eq!(
            infer_rank_from_lineage("k__Bacteria|p__P|c__C|o__O|f__F|g__G|s__S|t__SGB42"),
            "strain"
        );
    }

    #[test]
    fn writes_cami_header_and_rows() {
        let entries = vec![NormalizedEntry {
            identifier: "t__SGB42".to_owned(),
            lineage: "k__Bacteria|p__P|c__C|o__O|f__F|g__G|s__S|t__SGB42".to_owned(),
            abundance: 12.5,
            vertical_coverage: None,
            metadata: {
                let mut map = BTreeMap::new();
                map.insert(
                    "taxpath_ids".to_owned(),
                    "2|1224|1|2|3|4|5|987654".to_owned(),
                );
                map
            },
            source_format: "cami".to_owned(),
            source_tool: "metaphlan".to_owned(),
            source_version: "X".to_owned(),
            source_taxonomy: "GTDB".to_owned(),
        }];
        let mut out = Vec::new();
        write_cami(&mut out, "sample.profile", "X", &entries);
        let text = String::from_utf8(out).expect("invalid utf8");
        assert!(text.contains("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"));
        assert!(text.contains("987654\tstrain\t2|1224|1|2|3|4|5|987654\t"));
    }
}
