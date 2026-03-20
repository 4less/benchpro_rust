use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

/// NCBI taxonomy helper backed by `nodes.dmp` and `names.dmp`.
#[derive(Debug, Default)]
pub struct NcbiTaxdump {
    parent_by_taxid: HashMap<u64, u64>,
    rank_by_taxid: HashMap<u64, String>,
    scientific_name_by_taxid: HashMap<u64, String>,
}

impl NcbiTaxdump {
    /// Loads taxonomy from the default Benchpro cache directory.
    ///
    /// The cache root is `~/.benchpro/ncbi_taxdump` by default.
    /// Set `BENCHPRO_NCBI_DIR` to override this location.
    ///
    /// # Errors
    ///
    /// Returns an error if taxonomy files cannot be downloaded, extracted, or parsed.
    pub fn load_or_prepare() -> Result<Self, String> {
        let root = Self::taxdump_root()?;
        Self::load_or_prepare_from_root(&root)
    }

    /// Loads taxonomy from a specific directory.
    ///
    /// # Errors
    ///
    /// Returns an error if taxonomy files cannot be read or parsed.
    pub fn load_or_prepare_from_root(root: &Path) -> Result<Self, String> {
        fs::create_dir_all(root)
            .map_err(|err| format!("Cannot create taxdump directory '{}': {err}", root.display()))?;

        let nodes = root.join("nodes.dmp");
        let names = root.join("names.dmp");
        if !nodes.exists() || !names.exists() {
            Self::download_and_extract(root)?;
        }

        Self::from_dump_files(&nodes, &names)
    }

    /// Builds an NCBI lineage string from a taxid.
    ///
    /// # Returns
    ///
    /// A `|`-delimited lineage string with known scientific names.
    pub fn lineage_string_by_taxid(&self, taxid: u64) -> Option<String> {
        let mut lineage_ids = self.lineage_ids(taxid)?;
        lineage_ids.reverse();

        let names = lineage_ids
            .iter()
            .filter_map(|id| self.scientific_name_by_taxid.get(id))
            .map(|name| name.replace('|', " "))
            .collect::<Vec<_>>();

        if names.is_empty() {
            return None;
        }
        Some(names.join("|"))
    }

    /// Builds an NCBI lineage string restricted to:
    /// superkingdom, phylum, class, order, family, genus, species.
    ///
    /// Some taxonomies use `kingdom`/`domain` instead of `superkingdom`;
    /// those are accepted as fallback for the first rank.
    pub fn lineage_string_standard_ranks(&self, taxid: u64) -> Option<String> {
        let mut lineage_ids = self.lineage_ids(taxid)?;
        lineage_ids.reverse();

        let mut names = Vec::new();
        let rank_groups: [&[&str]; 7] = [
            &["superkingdom", "kingdom", "domain"],
            &["phylum"],
            &["class"],
            &["order"],
            &["family"],
            &["genus"],
            &["species"],
        ];

        for ranks in rank_groups {
            let id = lineage_ids.iter().find(|id| {
                self.rank_by_taxid
                    .get(id)
                    .is_some_and(|r| ranks.iter().any(|wanted| r == wanted))
            });
            if let Some(id) = id {
                if let Some(name) = self.scientific_name_by_taxid.get(id) {
                    names.push(name.replace('|', " "));
                }
            }
        }

        if names.is_empty() {
            return self.lineage_string_by_taxid(taxid);
        }
        Some(names.join("|"))
    }

    fn lineage_ids(&self, taxid: u64) -> Option<Vec<u64>> {
        if !self.parent_by_taxid.contains_key(&taxid) {
            return None;
        }

        let mut result = Vec::new();
        let mut current = taxid;
        let mut guard = 0usize;
        while guard < 512 {
            result.push(current);
            let parent = *self.parent_by_taxid.get(&current)?;
            if parent == current {
                break;
            }
            current = parent;
            guard += 1;
        }
        Some(result)
    }

    fn from_dump_files(nodes_path: &Path, names_path: &Path) -> Result<Self, String> {
        let mut parent_by_taxid = HashMap::new();
        let mut rank_by_taxid = HashMap::new();
        let mut scientific_name_by_taxid = HashMap::new();

        let nodes_file = File::open(nodes_path).map_err(|err| {
            format!("Cannot open nodes.dmp '{}': {err}", nodes_path.display())
        })?;
        for line in BufReader::new(nodes_file).lines() {
            let line = line.map_err(|err| {
                format!("Failed reading nodes.dmp '{}': {err}", nodes_path.display())
            })?;
            let fields = split_dmp_fields(&line);
            if fields.len() < 3 {
                continue;
            }
            let taxid = match fields[0].parse::<u64>() {
                Ok(value) => value,
                Err(_) => continue,
            };
            let parent = match fields[1].parse::<u64>() {
                Ok(value) => value,
                Err(_) => continue,
            };
            parent_by_taxid.insert(taxid, parent);
            rank_by_taxid.insert(taxid, fields[2].to_owned());
        }

        let names_file = File::open(names_path).map_err(|err| {
            format!("Cannot open names.dmp '{}': {err}", names_path.display())
        })?;
        for line in BufReader::new(names_file).lines() {
            let line = line.map_err(|err| {
                format!("Failed reading names.dmp '{}': {err}", names_path.display())
            })?;
            let fields = split_dmp_fields(&line);
            if fields.len() < 4 || fields[3] != "scientific name" {
                continue;
            }
            let taxid = match fields[0].parse::<u64>() {
                Ok(value) => value,
                Err(_) => continue,
            };
            scientific_name_by_taxid.insert(taxid, fields[1].to_owned());
        }

        Ok(Self {
            parent_by_taxid,
            rank_by_taxid,
            scientific_name_by_taxid,
        })
    }

    fn taxdump_root() -> Result<PathBuf, String> {
        if let Ok(path) = env::var("BENCHPRO_NCBI_DIR") {
            let path = path.trim();
            if !path.is_empty() {
                return Ok(PathBuf::from(path));
            }
        }

        let home = env::var("HOME")
            .map_err(|_| "HOME is not set and BENCHPRO_NCBI_DIR was not provided".to_owned())?;
        Ok(Path::new(&home).join(".benchpro").join("ncbi_taxdump"))
    }

    fn download_and_extract(root: &Path) -> Result<(), String> {
        let archive_path = root.join("new_taxdump.tar.gz");
        let url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz";

        let curl_status = Command::new("curl")
            .args(["-L", "-o"])
            .arg(&archive_path)
            .arg(url)
            .status()
            .map_err(|err| format!("Failed to run curl for NCBI taxdump download: {err}"))?;

        if !curl_status.success() {
            return Err(format!(
                "Downloading NCBI taxdump failed. Tried URL: {url}. \
                 You can manually download and extract nodes.dmp/names.dmp into '{}'.",
                root.display()
            ));
        }

        let tar_status = Command::new("tar")
            .args(["-xzf"])
            .arg(&archive_path)
            .args(["-C"])
            .arg(root)
            .args(["nodes.dmp", "names.dmp"])
            .status()
            .map_err(|err| format!("Failed to run tar for NCBI taxdump extraction: {err}"))?;

        if !tar_status.success() {
            return Err(format!(
                "Extracting NCBI taxdump failed. \
                 Ensure tar is installed and '{}' is a valid archive.",
                archive_path.display()
            ));
        }

        Ok(())
    }
}

fn split_dmp_fields(line: &str) -> Vec<&str> {
    line.split('|').map(str::trim).collect::<Vec<_>>()
}
