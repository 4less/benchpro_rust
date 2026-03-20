use regex::Regex;

/// Top-level profile format category used by `normalize --detect`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProfileFormatKind {
    Cami,
    Tool,
    Unknown,
}

impl ProfileFormatKind {
    /// Returns the machine-readable format string.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Cami => "cami",
            Self::Tool => "tool",
            Self::Unknown => "unknown",
        }
    }
}

/// Result returned by profile detection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DetectionResult {
    /// High-level format (`cami` or `tool`).
    pub format: ProfileFormatKind,
    /// Detected tool name (for example `motus`, `metaphlan`).
    pub tool: String,
    /// Detected tool or profile version, or `X` if unavailable.
    pub version: String,
}

impl DetectionResult {
    /// Creates a detection result with unknown tool and version.
    pub fn unknown() -> Self {
        Self {
            format: ProfileFormatKind::Unknown,
            tool: "X".to_owned(),
            version: "X".to_owned(),
        }
    }
}

/// Detects profile format and tool using content-based heuristics.
///
/// # Arguments
///
/// * `content` - Full profile file content
///
/// # Returns
///
/// Detection result containing format, tool, and version.
pub fn detect_profile(content: &str) -> DetectionResult {
    let lines = content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .collect::<Vec<_>>();

    if lines.is_empty() {
        return DetectionResult::unknown();
    }

    if is_cami(&lines) {
        return detect_cami(&lines);
    }

    detect_tool(&lines)
}

fn is_cami(lines: &[&str]) -> bool {
    lines
        .iter()
        .any(|line| line.starts_with("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"))
}

fn detect_cami(lines: &[&str]) -> DetectionResult {
    let joined = lines.join("\n").to_lowercase();
    let tool = if joined.contains("motus version") || joined.contains("motus profile") {
        "motus"
    } else if joined.contains("metaphlan") || joined.contains("#mpa_") {
        "metaphlan"
    } else {
        "X"
    };

    // Prefer explicit tool version when available, then CAMI @Version.
    let version = extract_motus_version(lines)
        .or_else(|| extract_generic_version(lines))
        .unwrap_or_else(|| "X".to_owned());

    DetectionResult {
        format: ProfileFormatKind::Cami,
        tool: tool.to_owned(),
        version,
    }
}

fn detect_tool(lines: &[&str]) -> DetectionResult {
    let first_non_comment = lines
        .iter()
        .copied()
        .find(|line| !line.starts_with('#'))
        .unwrap_or(lines[0]);

    if is_bracken(lines) {
        return DetectionResult {
            format: ProfileFormatKind::Tool,
            tool: "bracken".to_owned(),
            version: extract_keyword_version(lines, "bracken").unwrap_or_else(|| "X".to_owned()),
        };
    }

    if is_sylph(lines) {
        return DetectionResult {
            format: ProfileFormatKind::Tool,
            tool: "sylph".to_owned(),
            version: extract_keyword_version(lines, "sylph").unwrap_or_else(|| "X".to_owned()),
        };
    }

    if is_metaphlan_tool(lines) {
        return DetectionResult {
            format: ProfileFormatKind::Tool,
            tool: "metaphlan".to_owned(),
            version: extract_metaphlan_version(lines).unwrap_or_else(|| "X".to_owned()),
        };
    }

    if is_motus_relab(lines) {
        return DetectionResult {
            format: ProfileFormatKind::Tool,
            tool: "motus".to_owned(),
            version: extract_motus_tool_version(lines).unwrap_or_else(|| "X".to_owned()),
        };
    }

    if is_protal(lines) {
        return DetectionResult {
            format: ProfileFormatKind::Tool,
            tool: "protal".to_owned(),
            version: "X".to_owned(),
        };
    }

    if looks_like_two_col_lineage_abundance(first_non_comment) {
        return DetectionResult {
            format: ProfileFormatKind::Tool,
            tool: "mg-tk".to_owned(),
            version: "X".to_owned(),
        };
    }

    DetectionResult {
        format: ProfileFormatKind::Unknown,
        tool: "X".to_owned(),
        version: "X".to_owned(),
    }
}

fn is_bracken(lines: &[&str]) -> bool {
    lines.iter().any(|line| {
        let lower = line.to_ascii_lowercase();
        lower.starts_with("name\ttaxonomy_id\ttaxonomy_lvl")
            && lower.contains("fraction_total_reads")
    })
}

fn is_sylph(lines: &[&str]) -> bool {
    lines.iter().any(|line| {
        let lower = line.to_ascii_lowercase();
        lower.starts_with("clade_name\trelative_abundance\tsequence_abundance")
            && lower.contains("coverage")
    })
}

fn is_metaphlan_tool(lines: &[&str]) -> bool {
    lines.iter().any(|line| {
        let lower = line.to_ascii_lowercase();
        lower.starts_with("#clade_name\t") && lower.contains("relative_abundance")
            || lower == "#clade_name\trelative_abundance"
    })
}

fn is_protal(lines: &[&str]) -> bool {
    lines.iter().filter(|line| !line.starts_with('#')).take(8).any(|line| {
        let tokens = line.split('\t').collect::<Vec<_>>();
        if tokens.len() < 3 {
            return false;
        }
        let first_col = tokens[0].trim().to_ascii_lowercase();
        if first_col == "motu" || first_col.starts_with("motuv") {
            return false;
        }
        let lineage = tokens[1];
        let abundance_ok = tokens[2].trim().parse::<f64>().is_ok();
        abundance_ok
            && (lineage.contains("d__") || lineage.contains("k__"))
            && (lineage.contains('|') || lineage.contains(';'))
    })
}

fn is_motus_relab(lines: &[&str]) -> bool {
    lines.iter().any(|line| {
        let lower = line.to_ascii_lowercase();
        lower.starts_with("motu\ttaxonomy\t")
    })
}

fn looks_like_two_col_lineage_abundance(line: &str) -> bool {
    let tokens = line.split('\t').collect::<Vec<_>>();
    if tokens.len() < 2 {
        return false;
    }

    let lineage = tokens[0];
    let abundance_ok = tokens[1].parse::<f64>().is_ok();

    abundance_ok
        && (lineage.contains("d__") || lineage.contains("k__"))
        && (lineage.contains('|') || lineage.contains(';'))
}

fn extract_generic_version(lines: &[&str]) -> Option<String> {
    lines
        .iter()
        .find_map(|line| line.strip_prefix("@Version:"))
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .map(ToOwned::to_owned)
}

fn extract_motus_version(lines: &[&str]) -> Option<String> {
    let re = Regex::new(r"(?i)motus\s+version\s+([0-9]+(?:\.[0-9]+){0,3})").ok()?;
    lines.iter().find_map(|line| {
        re.captures(line)
            .and_then(|caps| caps.get(1).map(|m| m.as_str().to_owned()))
    })
}

fn extract_motus_tool_version(lines: &[&str]) -> Option<String> {
    let re = Regex::new(r"(?i)tool_version=([0-9]+(?:\.[0-9]+){0,3})").ok()?;
    lines.iter().find_map(|line| {
        re.captures(line)
            .and_then(|caps| caps.get(1).map(|m| m.as_str().to_owned()))
    })
}

fn extract_keyword_version(lines: &[&str], keyword: &str) -> Option<String> {
    let pattern = format!(
        r"(?i){}[^0-9]*([0-9]+(?:\.[0-9]+){{0,3}})",
        regex::escape(keyword)
    );
    let re = Regex::new(&pattern).ok()?;
    lines.iter().find_map(|line| {
        re.captures(line)
            .and_then(|caps| caps.get(1).map(|m| m.as_str().to_owned()))
    })
}

fn extract_metaphlan_version(lines: &[&str]) -> Option<String> {
    extract_keyword_version(lines, "metaphlan")
}

#[cfg(test)]
mod tests {
    use super::{detect_profile, ProfileFormatKind};

    #[test]
    fn detects_bracken_tool() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/bracken.profile"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Tool);
        assert_eq!(result.tool, "bracken");
    }

    #[test]
    fn detects_metaphlan_cami_and_version() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/metaphlan3013_camioutput.profile"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Cami);
        assert_eq!(result.tool, "metaphlan");
        assert_eq!(result.version, "0.10.0");
    }

    #[test]
    fn detects_motus_cami_and_prefers_tool_version() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/motus301_parenthesis.profile"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Cami);
        assert_eq!(result.tool, "motus");
        assert_eq!(result.version, "3.0.1");
    }

    #[test]
    fn detects_metaphlan_tool_format() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/metaphlan402.profile"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Tool);
        assert_eq!(result.tool, "metaphlan");
    }

    #[test]
    fn detects_mg_tk_like_two_column_profiles() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/test_data/profiles/predictions/MG-TK/GIT.sample.0.profile"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Tool);
        assert_eq!(result.tool, "mg-tk");
    }

    #[test]
    fn detects_sylph_profile() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/test_data/profiles/predictions/sylph/Gastrointestinal_1.profile"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Tool);
        assert_eq!(result.tool, "sylph");
    }

    #[test]
    fn detects_sylph_profile_example_repo_file() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/sylph090.profile"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Tool);
        assert_eq!(result.tool, "sylph");
    }

    #[test]
    fn detects_protal_profile() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/protal.profile"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Tool);
        assert_eq!(result.tool, "protal");
    }

    #[test]
    fn detects_motus404_relab_profile() {
        let content = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/profile_examples/motus404.relab"
        ));
        let result = detect_profile(content);
        assert_eq!(result.format, ProfileFormatKind::Tool);
        assert_eq!(result.tool, "motus");
        assert_eq!(result.version, "4.0.4");
    }
}
