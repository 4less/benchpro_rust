use std::collections::BTreeMap;

use crate::ncbi_taxdump::NcbiTaxdump;
use crate::normalize_detect::{DetectionResult, ProfileFormatKind};

/// One normalized record written by the `normalize` command.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedEntry {
    /// Taxon identifier, typically the right-most lineage token.
    pub identifier: String,
    /// Original lineage string when available.
    pub lineage: String,
    /// Reported abundance value.
    pub abundance: f64,
    /// Optional vertical coverage value.
    pub vertical_coverage: Option<f64>,
    /// Optional extensible metadata map for source-specific fields (for example tax IDs).
    pub metadata: BTreeMap<String, String>,
    /// Detected top-level source format (`cami`, `tool`, `unknown`).
    pub source_format: String,
    /// Detected source tool name.
    pub source_tool: String,
    /// Detected source version.
    pub source_version: String,
    /// Inferred source taxonomy label.
    pub source_taxonomy: String,
}

/// Loads and normalizes profile rows into a unified table shape.
///
/// # Arguments
///
/// * `content` - Full profile file content
/// * `detection` - Detection result from `detect_profile`
///
/// # Returns
///
/// Normalized entries with `identifier`, `lineage`, `abundance`, and optional `vertical_coverage`.
///
/// # Errors
///
/// Returns an error when the content cannot be parsed for the detected format.
pub fn load_normalized(
    content: &str,
    detection: &DetectionResult,
) -> Result<Vec<NormalizedEntry>, String> {
    load_normalized_for_kind(
        content,
        &detection.format,
        &detection.tool,
        &detection.version,
    )
}

/// Loads and normalizes profile rows from an explicit format/tool pair.
pub fn load_normalized_for_kind(
    content: &str,
    format: &ProfileFormatKind,
    tool: &str,
    version: &str,
) -> Result<Vec<NormalizedEntry>, String> {
    let mut entries = match format {
        ProfileFormatKind::Cami => load_cami(content),
        ProfileFormatKind::Tool => match tool {
            "bracken" => load_bracken(content),
            "sylph" => load_sylph(content),
            "metaphlan" => load_metaphlan_tool(content),
            "motus" => load_motus_relab(content),
            "protal" => load_protal_three_col(content),
            "mg-tk" => load_mg_tk_two_col(content),
            _ => load_generic_tool(content),
        },
        ProfileFormatKind::Unknown => load_generic_tool(content),
    }?;

    let taxonomy = infer_taxonomy(tool, &entries);
    for entry in &mut entries {
        entry.source_format = format.as_str().to_owned();
        entry.source_tool = tool.to_owned();
        entry.source_version = version.to_owned();
        entry.source_taxonomy = taxonomy.to_owned();
    }

    Ok(entries)
}

fn load_cami(content: &str) -> Result<Vec<NormalizedEntry>, String> {
    let lines = content.lines().collect::<Vec<_>>();
    let (header_idx, header) = lines
        .iter()
        .enumerate()
        .find_map(|(idx, line)| line.strip_prefix("@@").map(|h| (idx, h)))
        .ok_or_else(|| "CAMI header line starting with @@ not found".to_owned())?;

    let headers = header.split('\t').collect::<Vec<_>>();
    let taxpath_idx = headers.iter().position(|h| *h == "TAXPATH");
    let lineage_idx = headers
        .iter()
        .position(|h| *h == "TAXPATHSN")
        .ok_or_else(|| "CAMI TAXPATHSN column not found".to_owned())?;
    let abundance_idx = headers
        .iter()
        .position(|h| *h == "PERCENTAGE")
        .ok_or_else(|| "CAMI PERCENTAGE column not found".to_owned())?;

    let mut entries = Vec::with_capacity(lines.len().saturating_sub(header_idx + 1));
    for line in lines.iter().skip(header_idx + 1) {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with('@') {
            continue;
        }

        let tokens = trimmed.split('\t').collect::<Vec<_>>();
        if tokens.len() <= abundance_idx || tokens.len() <= lineage_idx {
            continue;
        }

        let lineage = trim_lineage_to_standard_ranks(tokens[lineage_idx].trim());
        let taxpath_ids = taxpath_idx.and_then(|idx| {
            tokens
                .get(idx)
                .map(|value| trim_taxpath_ids_to_match_lineage(value.trim(), &lineage))
        });
        let abundance = match tokens[abundance_idx].trim().parse::<f64>() {
            Ok(value) => value,
            Err(_) => continue,
        };

        entries.push(NormalizedEntry {
            identifier: rightmost_taxon(&lineage),
            lineage,
            abundance,
            vertical_coverage: None,
            metadata: metadata_with_taxpath_ids(taxpath_ids),
            source_format: String::new(),
            source_tool: String::new(),
            source_version: String::new(),
            source_taxonomy: String::new(),
        });
    }

    Ok(entries)
}

fn load_metaphlan_tool(content: &str) -> Result<Vec<NormalizedEntry>, String> {
    let lines = content.lines().collect::<Vec<_>>();
    let (header_idx, headers) = lines
        .iter()
        .enumerate()
        .find_map(|(idx, line)| {
            let cleaned = line.trim_start_matches('#').trim();
            if cleaned.to_ascii_lowercase().starts_with("clade_name\t")
                && cleaned.to_ascii_lowercase().contains("relative_abundance")
            {
                Some((idx, cleaned.split('\t').collect::<Vec<_>>()))
            } else {
                None
            }
        })
        .ok_or_else(|| {
            "MetaPhlAn header with clade_name and relative_abundance not found".to_owned()
        })?;

    let clade_idx = headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case("clade_name"))
        .ok_or_else(|| "MetaPhlAn clade_name column not found".to_owned())?;
    let abundance_idx = headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case("relative_abundance"))
        .ok_or_else(|| "MetaPhlAn relative_abundance column not found".to_owned())?;
    let mut entries = Vec::with_capacity(lines.len().saturating_sub(header_idx + 1));
    for line in lines.iter().skip(header_idx + 1) {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let tokens = trimmed.split('\t').collect::<Vec<_>>();
        if tokens.len() <= clade_idx || tokens.len() <= abundance_idx {
            continue;
        }

        let lineage = trim_lineage_to_standard_ranks(tokens[clade_idx].trim());
        let abundance = match tokens[abundance_idx].trim().parse::<f64>() {
            Ok(value) => value,
            Err(_) => continue,
        };

        entries.push(NormalizedEntry {
            identifier: rightmost_taxon(&lineage),
            lineage,
            abundance,
            vertical_coverage: None,
            metadata: BTreeMap::new(),
            source_format: String::new(),
            source_tool: String::new(),
            source_version: String::new(),
            source_taxonomy: String::new(),
        });
    }

    Ok(entries)
}

fn load_sylph(content: &str) -> Result<Vec<NormalizedEntry>, String> {
    let lines = content.lines().collect::<Vec<_>>();
    let (header_idx, headers) = lines
        .iter()
        .enumerate()
        .find_map(|(idx, line)| {
            let trimmed = line.trim();
            if trimmed
                .to_ascii_lowercase()
                .starts_with("clade_name\trelative_abundance\tsequence_abundance")
            {
                Some((idx, trimmed.split('\t').collect::<Vec<_>>()))
            } else {
                None
            }
        })
        .ok_or_else(|| "Sylph header not found".to_owned())?;

    let clade_idx = headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case("clade_name"))
        .ok_or_else(|| "Sylph clade_name column not found".to_owned())?;
    let abundance_idx = headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case("relative_abundance"))
        .ok_or_else(|| "Sylph relative_abundance column not found".to_owned())?;
    let coverage_idx = headers.iter().position(|h| {
        h.to_ascii_lowercase().contains("coverage")
            && h.to_ascii_lowercase().contains("strain-level")
    });

    let mut entries = Vec::with_capacity(lines.len().saturating_sub(header_idx + 1));
    for line in lines.iter().skip(header_idx + 1) {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let tokens = trimmed.split('\t').collect::<Vec<_>>();
        if tokens.len() <= clade_idx || tokens.len() <= abundance_idx {
            continue;
        }

        let lineage = trim_lineage_to_standard_ranks(tokens[clade_idx].trim());
        let abundance = match tokens[abundance_idx].trim().parse::<f64>() {
            Ok(value) => value,
            Err(_) => continue,
        };

        let vertical_coverage = coverage_idx.and_then(|idx| {
            tokens.get(idx).and_then(|value| {
                let normalized = value.trim();
                if normalized.eq_ignore_ascii_case("na") || normalized.is_empty() {
                    None
                } else {
                    normalized.parse::<f64>().ok()
                }
            })
        });

        entries.push(NormalizedEntry {
            identifier: rightmost_taxon(&lineage),
            lineage,
            abundance,
            vertical_coverage,
            metadata: BTreeMap::new(),
            source_format: String::new(),
            source_tool: String::new(),
            source_version: String::new(),
            source_taxonomy: String::new(),
        });
    }

    Ok(entries)
}

fn load_bracken(content: &str) -> Result<Vec<NormalizedEntry>, String> {
    let mut lines = content.lines();
    let header = lines
        .next()
        .ok_or_else(|| "Bracken file is empty".to_owned())?
        .split('\t')
        .collect::<Vec<_>>();

    let name_idx = header
        .iter()
        .position(|h| h.eq_ignore_ascii_case("name"))
        .ok_or_else(|| "Bracken name column not found".to_owned())?;
    let taxid_idx = header
        .iter()
        .position(|h| h.eq_ignore_ascii_case("taxonomy_id"))
        .ok_or_else(|| "Bracken taxonomy_id column not found".to_owned())?;
    let abundance_idx = header
        .iter()
        .position(|h| h.eq_ignore_ascii_case("fraction_total_reads"))
        .ok_or_else(|| "Bracken fraction_total_reads column not found".to_owned())?;

    // Bracken reports taxids but not full lineage strings, so we resolve lineages from NCBI taxdump.
    let ncbi = NcbiTaxdump::load_or_prepare()?;

    let mut entries = Vec::new();
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let tokens = trimmed.split('\t').collect::<Vec<_>>();
        if tokens.len() <= name_idx || tokens.len() <= abundance_idx || tokens.len() <= taxid_idx {
            continue;
        }

        let name = tokens[name_idx].trim();
        let taxid = tokens[taxid_idx].trim().parse::<u64>().ok();
        let lineage = taxid
            .and_then(|value| ncbi.lineage_string_standard_ranks(value))
            .unwrap_or_else(|| name.to_owned());
        let abundance = match tokens[abundance_idx].trim().parse::<f64>() {
            Ok(value) => value,
            Err(_) => continue,
        };

        entries.push(NormalizedEntry {
            identifier: rightmost_taxon(&lineage),
            lineage,
            abundance,
            vertical_coverage: None,
            metadata: metadata_with_taxpath_ids(taxid.map(|value| value.to_string())),
            source_format: String::new(),
            source_tool: String::new(),
            source_version: String::new(),
            source_taxonomy: String::new(),
        });
    }

    Ok(entries)
}

fn load_mg_tk_two_col(content: &str) -> Result<Vec<NormalizedEntry>, String> {
    let mut entries = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let tokens = trimmed.split('\t').collect::<Vec<_>>();
        if tokens.len() < 2 {
            continue;
        }

        let lineage = trim_lineage_to_standard_ranks(tokens[0].trim());
        let abundance = match tokens[1].trim().parse::<f64>() {
            Ok(value) => value,
            Err(_) => continue,
        };

        entries.push(NormalizedEntry {
            identifier: rightmost_taxon(&lineage),
            lineage,
            abundance,
            vertical_coverage: None,
            metadata: BTreeMap::new(),
            source_format: String::new(),
            source_tool: String::new(),
            source_version: String::new(),
            source_taxonomy: String::new(),
        });
    }

    if entries.is_empty() {
        return Err("Unable to parse MG-TK style two-column profile".to_owned());
    }

    Ok(entries)
}

fn load_protal_three_col(content: &str) -> Result<Vec<NormalizedEntry>, String> {
    let mut entries = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let tokens = trimmed.split('\t').collect::<Vec<_>>();
        if tokens.len() < 3 {
            continue;
        }

        let lineage = trim_lineage_to_standard_ranks(tokens[1].trim());
        let abundance = match tokens[2].trim().parse::<f64>() {
            Ok(value) => value,
            Err(_) => continue,
        };

        entries.push(NormalizedEntry {
            identifier: rightmost_taxon(&lineage),
            lineage,
            abundance,
            vertical_coverage: None,
            metadata: BTreeMap::new(),
            source_format: String::new(),
            source_tool: String::new(),
            source_version: String::new(),
            source_taxonomy: String::new(),
        });
    }

    if entries.is_empty() {
        return Err("Unable to parse Protal style three-column profile".to_owned());
    }

    Ok(entries)
}

fn load_motus_relab(content: &str) -> Result<Vec<NormalizedEntry>, String> {
    let mut entries = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let tokens = trimmed.split('\t').collect::<Vec<_>>();
        if tokens.len() < 3 {
            continue;
        }
        if tokens[0].eq_ignore_ascii_case("mOTU") && tokens[1].eq_ignore_ascii_case("taxonomy") {
            continue;
        }

        let lineage = trim_lineage_to_standard_ranks(tokens[1].trim());
        let abundance = tokens[2..]
            .iter()
            .find_map(|value| value.trim().parse::<f64>().ok());
        let Some(abundance) = abundance else {
            continue;
        };

        entries.push(NormalizedEntry {
            identifier: rightmost_taxon(&lineage),
            lineage,
            abundance,
            vertical_coverage: None,
            metadata: BTreeMap::new(),
            source_format: String::new(),
            source_tool: String::new(),
            source_version: String::new(),
            source_taxonomy: String::new(),
        });
    }

    if entries.is_empty() {
        return Err("Unable to parse mOTUs relab profile".to_owned());
    }

    Ok(entries)
}

fn load_generic_tool(content: &str) -> Result<Vec<NormalizedEntry>, String> {
    if let Ok(entries) = load_metaphlan_tool(content) {
        if !entries.is_empty() {
            return Ok(entries);
        }
    }

    if let Ok(entries) = load_mg_tk_two_col(content) {
        if !entries.is_empty() {
            return Ok(entries);
        }
    }

    Err("Could not parse profile with generic tool loaders".to_owned())
}

fn rightmost_taxon(lineage: &str) -> String {
    let delimiter = if lineage.contains('|') {
        '|'
    } else if lineage.contains(';') {
        ';'
    } else {
        return lineage.trim().to_owned();
    };

    lineage
        .split(delimiter)
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .last()
        .unwrap_or("")
        .to_owned()
}

fn trim_lineage_to_standard_ranks(lineage: &str) -> String {
    let delimiter = if lineage.contains('|') {
        '|'
    } else if lineage.contains(';') {
        ';'
    } else {
        return lineage.trim().to_owned();
    };

    let tokens = lineage
        .split(delimiter)
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .collect::<Vec<_>>();

    if tokens.is_empty() {
        return String::new();
    }

    // Prefix-aware lineages (for example d__/k__, p__, c__, ... plus optional t__).
    let mut by_rank: [Option<&str>; 8] = [None, None, None, None, None, None, None, None];
    let mut saw_prefix = false;
    for token in &tokens {
        if let Some((prefix, _)) = token.split_once("__") {
            saw_prefix = true;
            let index = match prefix {
                "d" | "k" | "sk" => Some(0),
                "p" => Some(1),
                "c" => Some(2),
                "o" => Some(3),
                "f" => Some(4),
                "g" => Some(5),
                "s" => Some(6),
                "t" => Some(7),
                _ => None,
            };
            if let Some(idx) = index {
                by_rank[idx] = Some(token);
            }
        }
    }

    if saw_prefix {
        let kept = by_rank.into_iter().flatten().collect::<Vec<_>>();
        if kept.is_empty() {
            return tokens.join(&delimiter.to_string());
        }
        return kept.join(&delimiter.to_string());
    }

    // Fallback for unprefixed lineage strings: keep at most eight ranks
    // (superkingdom..species plus optional strain).
    if tokens.len() >= 9 {
        return tokens[..8].join(&delimiter.to_string());
    }

    tokens.join(&delimiter.to_string())
}

fn trim_taxpath_ids_to_match_lineage(taxpath_ids: &str, lineage: &str) -> String {
    let lineage_len = lineage
        .split(['|', ';'])
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .count();
    if lineage_len == 0 {
        return String::new();
    }

    let tokens = taxpath_ids
        .split('|')
        .map(str::trim)
        .filter(|token| !token.is_empty())
        .collect::<Vec<_>>();
    if tokens.is_empty() {
        return String::new();
    }

    let keep = tokens.len().min(lineage_len);
    tokens[..keep].join("|")
}

fn metadata_with_taxpath_ids(taxpath_ids: Option<String>) -> BTreeMap<String, String> {
    let mut metadata = BTreeMap::new();
    if let Some(value) = taxpath_ids.filter(|value| !value.is_empty()) {
        metadata.insert("taxpath_ids".to_owned(), value);
    }
    metadata
}

fn infer_taxonomy(tool: &str, entries: &[NormalizedEntry]) -> &'static str {
    if tool.eq_ignore_ascii_case("bracken") {
        return "NCBI";
    }
    if tool.eq_ignore_ascii_case("metaphlan") {
        return "ChocoPhlAn";
    }
    if entries.iter().any(|entry| {
        entry.lineage.contains("d__")
            || entry.lineage.contains("k__")
            || entry.lineage.contains("p__")
    }) {
        return "GTDB";
    }
    if entries
        .iter()
        .any(|entry| entry.metadata.contains_key("taxpath_ids"))
    {
        return "NCBI";
    }
    "Custom"
}

#[cfg(test)]
mod tests {
    use crate::normalize_detect::detect_profile;

    use super::{load_normalized, trim_lineage_to_standard_ranks};

    #[test]
    fn loads_cami_content() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/metaphlan3013_camioutput.profile"
        ));
        let detection = detect_profile(content);
        let entries = load_normalized(content, &detection).expect("failed to load cami");
        assert!(!entries.is_empty());
        assert!(entries[0].lineage.contains("Bacteria"));
    }

    #[test]
    fn loads_sylph_with_vertical_coverage_when_present() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/test_data/profiles/predictions/sylph/Gastrointestinal_1.profile"
        ));
        let detection = detect_profile(content);
        let entries = load_normalized(content, &detection).expect("failed to load sylph");
        assert!(!entries.is_empty());
        assert!(entries
            .iter()
            .any(|entry| entry.vertical_coverage.is_some()));
        assert!(entries.iter().any(|entry| entry.lineage.contains("t__")));
    }

    #[test]
    fn loads_sylph_profile_example_repo_file() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/sylph090.profile"
        ));
        let detection = detect_profile(content);
        let entries = load_normalized(content, &detection).expect("failed to load sylph example");
        assert!(!entries.is_empty());
        assert!(entries
            .iter()
            .any(|entry| entry.vertical_coverage.is_some()));
    }

    #[test]
    fn keeps_metaphlan_t_rank_when_present() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/metaphlan402.profile"
        ));
        let detection = detect_profile(content);
        let entries = load_normalized(content, &detection).expect("failed to load metaphlan");
        assert!(entries.iter().any(|entry| entry.lineage.contains("t__")));
        assert!(entries
            .iter()
            .all(|entry| !entry.metadata.contains_key("taxpath_ids")));
    }

    #[test]
    fn trims_to_species_plus_optional_t_rank() {
        let lineage = "k__Bacteria|p__P|c__C|o__O|f__F|g__G|s__S|t__SGB1|x__ignored";
        let trimmed = trim_lineage_to_standard_ranks(lineage);
        assert_eq!(trimmed, "k__Bacteria|p__P|c__C|o__O|f__F|g__G|s__S|t__SGB1");
    }

    #[test]
    fn loads_protal_three_column_profile() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/protal.profile"
        ));
        let detection = detect_profile(content);
        let entries = load_normalized(content, &detection).expect("failed to load protal");
        assert!(!entries.is_empty());
        assert!(entries[0].lineage.contains("d__Bacteria"));
    }

    #[test]
    fn loads_motus404_relab_profile() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/motus404.relab"
        ));
        let detection = detect_profile(content);
        let entries = load_normalized(content, &detection).expect("failed to load motus404 relab");
        assert!(!entries.is_empty());
        assert!(entries[0].lineage.contains("d__Bacteria"));
    }
}
