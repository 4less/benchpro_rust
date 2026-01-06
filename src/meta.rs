use calamine::{Reader, Xlsx};
use itertools::{izip, Itertools};
use polars::{
    error::{PolarsError, PolarsResult},
    frame::DataFrame,
    io::SerReader,
    prelude::{
        col, lit, AnyValue, AsString, CsvReadOptions, DataFrameJoinOps, Expr, IntoLazy, NamedFrom,
        PlSmallStr, SortOptions,
    },
    series::{Series, SeriesIter},
};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    str::FromStr,
};
use strum::{EnumIter, IntoEnumIterator};
use log::debug;
use regex::Regex;

use crate::utils::workbook_to_dataframe;

/// Errors produced while parsing or validating meta files.
#[derive(thiserror::Error, Debug, Clone)]
pub enum MetaError {
    #[error("{0}")]
    MissingColumns(String),
    #[error("{0}")]
    Empty(String),
    #[error("{0}")]
    DataError(String),
}

pub type MetaResult = Result<Meta, MetaError>;

#[derive(Default)]
/// Parsed meta table and expanded entries.
pub struct Meta {
    pub raw: DataFrame,
    pub entries: Vec<MetaEntry>,
}

/// Column name constants for supported meta formats.
pub struct MetaColumnStrings;

impl MetaColumnStrings {
    pub const ID: &str = "ID";
    pub const SAMPLE: &str = "Sample";
    pub const DATASET: &str = "Dataset";
    pub const TOOL: &str = "Tool";
    pub const TAXONOMY: &str = "Taxonomy";
    pub const PROFILE: &str = "Profile";
    pub const PROFILE_COLUMNS: &str = "ProfileColumns";
    pub const GOLDSTD: &str = "GoldStd";
    pub const GOLDSTD_COLUMNS: &str = "GoldStdColumns";
    pub const GOLDSTD_TREE: &str = "GoldStdTree";
    pub const AVAILABLE_TAXA: &str = "AvailableSpecies";
    pub const PROFILE_REGEX: &str = "ProfileRegex";
    pub const GOLDSTD_REGEX: &str = "GoldStdRegex";
}

/// Enumeration of supported meta columns.
#[derive(Debug, EnumIter, Clone)]
pub enum MetaColumn {
    ID,
    Sample,
    Dataset,
    Tool,
    Taxonomy,
    Profile,
    ProfileColumns,
    GoldStd,
    GoldStdColumns,
    GoldStdTree,
    AvailableTaxa,
}

impl MetaColumn {
    /// Returns the canonical column name string.
    ///
    /// # Returns
    ///
    /// Column name as a string slice.
    pub fn to_str(&self) -> &str {
        match self {
            MetaColumn::ID => MetaColumnStrings::ID,
            MetaColumn::Sample => MetaColumnStrings::SAMPLE,
            MetaColumn::Dataset => MetaColumnStrings::DATASET,
            MetaColumn::Tool => MetaColumnStrings::TOOL,
            MetaColumn::Taxonomy => MetaColumnStrings::TAXONOMY,
            MetaColumn::Profile => MetaColumnStrings::PROFILE,
            MetaColumn::ProfileColumns => MetaColumnStrings::PROFILE_COLUMNS,
            MetaColumn::GoldStd => MetaColumnStrings::GOLDSTD,
            MetaColumn::GoldStdColumns => MetaColumnStrings::GOLDSTD_COLUMNS,
            MetaColumn::GoldStdTree => MetaColumnStrings::GOLDSTD_TREE,
            MetaColumn::AvailableTaxa => MetaColumnStrings::AVAILABLE_TAXA,
        }
    }

    /// Parses a column enum from a column name.
    ///
    /// # Arguments
    ///
    /// * `str` - Column name as provided in the meta file
    ///
    /// # Returns
    ///
    /// `Some(MetaColumn)` for recognized columns, otherwise `None`.
    pub fn from_string(str: &str) -> Option<Self> {
        match str {
            MetaColumnStrings::ID => Some(Self::ID),
            MetaColumnStrings::SAMPLE => Some(Self::Sample),
            MetaColumnStrings::DATASET => Some(Self::Dataset),
            MetaColumnStrings::TOOL => Some(Self::Tool),
            MetaColumnStrings::TAXONOMY => Some(Self::Taxonomy),
            MetaColumnStrings::PROFILE => Some(Self::Profile),
            MetaColumnStrings::PROFILE_COLUMNS => Some(Self::ProfileColumns),
            MetaColumnStrings::GOLDSTD => Some(Self::GoldStd),
            MetaColumnStrings::GOLDSTD_COLUMNS => Some(Self::GoldStdColumns),
            MetaColumnStrings::GOLDSTD_TREE => Some(Self::GoldStdTree),
            MetaColumnStrings::AVAILABLE_TAXA => Some(Self::AvailableTaxa),
            "AvailableSpecies" => Some(Self::AvailableTaxa),
            _ => None,
        }
    }
}

impl From<MetaColumn> for String {
    fn from(value: MetaColumn) -> Self {
        value.to_str().to_string()
    }
}

impl From<MetaColumn> for PlSmallStr {
    fn from(value: MetaColumn) -> Self {
        PlSmallStr::from_str(value.to_str())
    }
}

impl Into<Expr> for MetaColumn {
    fn into(self) -> Expr {
        Expr::from(self.to_str())
    }
}

/// Single meta row expanded into a resolved input pair.
#[derive(Default, Clone, Debug)]
pub struct MetaEntry {
    pub id: String,
    pub sample: String,
    pub dataset: Option<String>,
    pub tool: String,
    pub taxonomy: String,
    pub profile: PathBuf,
    pub profile_columns: Option<String>,
    pub goldstd: PathBuf,
    pub goldstd_columns: Option<String>,
    pub goldstd_tree: Option<PathBuf>,
    pub taxa_list: Option<PathBuf>,
}

impl Meta {
    const REQUIRED_FIELDS: &[&str] = &["ID", "Sample", "Tool", "Profile", "GoldStd"];
    const TOTAL_FIELDS: &[&str] = &[
        "ID",
        "Sample",
        "Dataset",
        "Tool",
        "Taxonomy",
        "Profile",
        "GoldStd",
        "GoldStdColumns",
        "ProfileFormat",
    ];

    /// Builds a `Meta` from an in-memory DataFrame.
    ///
    /// # Arguments
    ///
    /// * `meta` - DataFrame with meta columns
    ///
    /// # Returns
    ///
    /// Parsed `Meta` with expanded entries.
    ///
    /// # Errors
    ///
    /// Returns `MetaError` when required columns are missing or data is invalid.
    pub fn from_polars_df(meta: DataFrame) -> MetaResult {
        if Self::is_matrix_format(&meta) {
            return Self::from_matrix_df(meta);
        }
        let mut entries = vec![MetaEntry::default(); meta.height()];

        debug!("Height: {} .. {}", meta.height(), entries.len());

        meta.get_column_names()
            .iter()
            .map(|x| x.to_string())
            .for_each(|name| {
                debug!("{} -> {:?}", name, MetaColumn::from_string(&name));
            });

        MetaColumn::iter().for_each(|col| {
            let column_series = meta.column(col.to_str());

            if let Ok(cs) = column_series {
                cs.str()
                    .unwrap()
                    .iter()
                    .enumerate()
                    .for_each(|(row_index, entry)| {
                        let entry = entry.expect("Expect string ");
                        let meta_entry = &mut entries[row_index];
                        match col {
                            MetaColumn::ID => meta_entry.id = entry.to_string(),
                            MetaColumn::Sample => meta_entry.sample = entry.to_string(),
                            MetaColumn::Dataset => meta_entry.dataset = Some(entry.to_string()),
                            MetaColumn::Tool => meta_entry.tool = entry.to_string(),
                            MetaColumn::Taxonomy => meta_entry.taxonomy = entry.to_string(),
                            MetaColumn::Profile => {
                                meta_entry.profile = PathBuf::from(entry.to_string())
                            }
                            MetaColumn::ProfileColumns => {
                                meta_entry.profile_columns = Some(entry.to_string())
                            }
                            MetaColumn::GoldStd => {
                                meta_entry.goldstd = PathBuf::from(entry.to_string())
                            }
                            MetaColumn::GoldStdColumns => {
                                meta_entry.goldstd_columns = Some(entry.to_string())
                            }
                            MetaColumn::GoldStdTree => {
                                meta_entry.goldstd_tree = Some(PathBuf::from(entry.to_string()))
                            }
                            MetaColumn::AvailableTaxa => {
                                meta_entry.taxa_list = match entry {
                                    entry if entry == "NA" || entry == "" => None,
                                    _ => Some(PathBuf::from(entry.to_string())),
                                }
                            }
                        }
                    });
            }
        });

        Ok(Self { raw: meta, entries })
    }

    pub fn left_join_to(
        &self,
        left: &DataFrame,
        columns: &[MetaColumn],
        remove_duplicates: bool,
    ) -> PolarsResult<DataFrame> {
        let mut columns = columns
            .iter()
            .map(|x| x.to_str())
            .collect::<HashSet<&str>>();
        columns.insert(MetaColumn::ID.to_str());

        if remove_duplicates {
            let present_columns = left
                .get_column_names()
                .iter()
                .map(|x| x.as_str())
                .collect::<HashSet<&str>>();
            columns = columns
                .difference(&present_columns)
                .map(|&x| x)
                .collect::<HashSet<_>>();
            columns.insert(MetaColumn::ID.to_str());
        }

        // let columns = columns.
        left.left_join(
            &self.raw.select(columns).expect("Cannot subset df"),
            [MetaColumn::ID.to_str()],
            [MetaColumn::ID.to_str()],
        )
    }

    fn df_from_text(path: impl AsRef<Path>) -> DataFrame {
        let df = CsvReadOptions::default()
            .try_into_reader_with_file_path(Some(path.as_ref().to_path_buf()))
            .expect("Cannot read file")
            .finish()
            .unwrap();
        df
    }

    fn df_from_xlsx(path: impl AsRef<Path>) -> DataFrame {
        // Open the Excel file
        let mut workbook: Xlsx<_> = calamine::open_workbook(path).expect("Cannot open Excel file");

        let df = workbook_to_dataframe(&mut workbook);

        df.unwrap()
    }

    /// Validates required columns and consistency constraints.
    ///
    /// # Returns
    ///
    /// `Meta` if validation succeeds.
    ///
    /// # Errors
    ///
    /// Returns `MetaError` when columns are missing or inconsistent.
    pub fn validate(self) -> MetaResult {
        type C = MetaColumn;
        // let colnames =
        let required_fields_present = Self::REQUIRED_FIELDS.iter().all(|&field| {
            self.raw
                .get_column_names()
                .into_iter()
                .map(|name| name.to_string())
                .collect::<Vec<_>>()
                .contains(&field.to_string())
        });
        debug!(
            "Columns({}) {:?}",
            required_fields_present,
            self.raw.get_column_names()
        );

        if !required_fields_present {
            return Err(MetaError::MissingColumns(
                "Not all columns are present".to_string(),
            ));
        }

        let check_goldstd_taxonomy = self
            .raw
            .clone()
            .lazy()
            .group_by([C::GoldStd])
            .agg([col(C::Taxonomy).n_unique().alias("UniqueCount")])
            .with_column(col("UniqueCount").eq(lit(1)).alias("IsUnique"))
            .collect()
            .expect("no error");
        let all_valid = check_goldstd_taxonomy
            .column("IsUnique")
            .unwrap()
            .bool()
            .unwrap()
            .all();

        if !all_valid {
            let error_df = check_goldstd_taxonomy
                .lazy()
                .filter(col("IsUnique").eq(lit(false)))
                .collect()
                .unwrap();
            let series = error_df
                .column(C::GoldStd.to_str())
                .unwrap()
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>();
            return Err(MetaError::DataError(format!("GoldStd profiles occurring more than once in meta must always have the same Taxonomy \n{:?}", series)));
        }

        Ok(self)
    }

    fn is_matrix_format(meta: &DataFrame) -> bool {
        let names = meta.get_column_names();
        names
            .iter()
            .any(|name| name.as_str() == MetaColumnStrings::PROFILE_REGEX)
            && names
                .iter()
                .any(|name| name.as_str() == MetaColumnStrings::GOLDSTD_REGEX)
    }

    fn read_matrix_headers(path: &Path) -> Result<Vec<String>, MetaError> {
        let file = File::open(path).map_err(|e| {
            MetaError::DataError(format!("Cannot open abundance matrix '{}': {}", path.display(), e))
        })?;
        let mut reader = BufReader::new(file);
        let mut header = String::new();
        reader.read_line(&mut header).map_err(|e| {
            MetaError::DataError(format!(
                "Cannot read header from abundance matrix '{}': {}",
                path.display(),
                e
            ))
        })?;
        let header = header.trim_end_matches(['\n', '\r']);
        if header.is_empty() {
            return Err(MetaError::DataError(format!(
                "Abundance matrix '{}' has empty header",
                path.display()
            )));
        }
        Ok(header.split('\t').map(|s| s.to_string()).collect())
    }

    fn match_headers_with_regex(
        headers: &[String],
        raw_regex: &str,
        label: &str,
    ) -> Result<std::collections::HashMap<String, String>, MetaError> {
        let raw_regex = raw_regex.trim();
        if raw_regex.is_empty() || raw_regex == "NA" {
            return Err(MetaError::DataError(format!(
                "{} is empty or NA in matrix meta",
                label
            )));
        }

        let candidates = if raw_regex.contains("\\\\") {
            vec![raw_regex.to_string(), raw_regex.replace("\\\\", "\\")]
        } else {
            vec![raw_regex.to_string()]
        };

        for candidate in candidates {
            let regex = Regex::new(&candidate).map_err(|e| {
                MetaError::DataError(format!(
                    "Invalid {} regex '{}': {}",
                    label, candidate, e
                ))
            })?;
            let mut matches = std::collections::HashMap::new();

            for header in headers.iter().skip(1) {
                if let Some(captures) = regex.captures(header) {
                    let capture = captures.get(1).ok_or_else(|| {
                        MetaError::DataError(format!(
                            "{} regex '{}' matched '{}' but has no capture group",
                            label, candidate, header
                        ))
                    })?;
                    let key = capture.as_str().to_string();
                    if matches.insert(key.clone(), header.clone()).is_some() {
                        return Err(MetaError::DataError(format!(
                            "{} regex '{}' produced duplicate capture key '{}'",
                            label, candidate, key
                        )));
                    }
                }
            }

            if !matches.is_empty() {
                return Ok(matches);
            }
        }

        Err(MetaError::DataError(format!(
            "No columns match {} regex '{}'",
            label, raw_regex
        )))
    }

    fn normalize_optional_str(value: Option<&str>) -> Option<String> {
        let value = value?.trim();
        if value.is_empty() || value == "NA" {
            None
        } else {
            Some(value.to_string())
        }
    }

    fn from_matrix_df(meta: DataFrame) -> MetaResult {
        let id_series = meta
            .column(MetaColumnStrings::ID)
            .map_err(|e| MetaError::MissingColumns(format!("Missing ID column: {}", e)))?
            .str()
            .map_err(|e| MetaError::DataError(format!("Invalid ID column: {}", e)))?;
        let dataset_series = meta
            .column(MetaColumnStrings::DATASET)
            .ok()
            .and_then(|s| s.str().ok());
        let tool_series = meta
            .column(MetaColumnStrings::TOOL)
            .map_err(|e| MetaError::MissingColumns(format!("Missing Tool column: {}", e)))?
            .str()
            .map_err(|e| MetaError::DataError(format!("Invalid Tool column: {}", e)))?;
        let taxonomy_series = meta
            .column(MetaColumnStrings::TAXONOMY)
            .map_err(|e| MetaError::MissingColumns(format!("Missing Taxonomy column: {}", e)))?
            .str()
            .map_err(|e| MetaError::DataError(format!("Invalid Taxonomy column: {}", e)))?;
        let profile_series = meta
            .column(MetaColumnStrings::PROFILE)
            .map_err(|e| MetaError::MissingColumns(format!("Missing Profile column: {}", e)))?
            .str()
            .map_err(|e| MetaError::DataError(format!("Invalid Profile column: {}", e)))?;
        let profile_regex_series = meta
            .column(MetaColumnStrings::PROFILE_REGEX)
            .map_err(|e| MetaError::MissingColumns(format!("Missing ProfileRegex column: {}", e)))?
            .str()
            .map_err(|e| MetaError::DataError(format!("Invalid ProfileRegex column: {}", e)))?;
        let goldstd_series = meta
            .column(MetaColumnStrings::GOLDSTD)
            .map_err(|e| MetaError::MissingColumns(format!("Missing GoldStd column: {}", e)))?
            .str()
            .map_err(|e| MetaError::DataError(format!("Invalid GoldStd column: {}", e)))?;
        let goldstd_regex_series = meta
            .column(MetaColumnStrings::GOLDSTD_REGEX)
            .map_err(|e| MetaError::MissingColumns(format!("Missing GoldStdRegex column: {}", e)))?
            .str()
            .map_err(|e| MetaError::DataError(format!("Invalid GoldStdRegex column: {}", e)))?;
        let goldstd_tree_series = meta
            .column(MetaColumnStrings::GOLDSTD_TREE)
            .ok()
            .and_then(|s| s.str().ok());
        let taxa_series = meta
            .column(MetaColumnStrings::AVAILABLE_TAXA)
            .ok()
            .and_then(|s| s.str().ok());

        let mut entries = Vec::new();

        let mut ids = Vec::new();
        let mut samples = Vec::new();
        let mut datasets = Vec::new();
        let mut tools = Vec::new();
        let mut taxonomies = Vec::new();
        let mut profiles = Vec::new();
        let mut goldstds = Vec::new();
        let mut goldstd_trees = Vec::new();
        let mut taxa_lists = Vec::new();

        let mut header_cache: std::collections::HashMap<String, Vec<String>> =
            std::collections::HashMap::new();

        for row_index in 0..meta.height() {
            let base_id = id_series.get(row_index).ok_or_else(|| {
                MetaError::DataError("Missing ID value in matrix meta".to_string())
            })?;
            let tool = tool_series.get(row_index).ok_or_else(|| {
                MetaError::DataError("Missing Tool value in matrix meta".to_string())
            })?;
            let taxonomy = taxonomy_series.get(row_index).ok_or_else(|| {
                MetaError::DataError("Missing Taxonomy value in matrix meta".to_string())
            })?;
            let profile_path = profile_series.get(row_index).ok_or_else(|| {
                MetaError::DataError("Missing Profile path in matrix meta".to_string())
            })?.trim();
            let profile_regex = profile_regex_series.get(row_index).ok_or_else(|| {
                MetaError::DataError("Missing ProfileRegex in matrix meta".to_string())
            })?;
            let goldstd_path = goldstd_series.get(row_index).ok_or_else(|| {
                MetaError::DataError("Missing GoldStd path in matrix meta".to_string())
            })?.trim();
            let goldstd_regex = goldstd_regex_series.get(row_index).ok_or_else(|| {
                MetaError::DataError("Missing GoldStdRegex in matrix meta".to_string())
            })?;

            let dataset = dataset_series
                .as_ref()
                .and_then(|s| s.get(row_index));
            let goldstd_tree = goldstd_tree_series
                .as_ref()
                .and_then(|s| s.get(row_index));
            let taxa_list = taxa_series.as_ref().and_then(|s| s.get(row_index));

            let profile_headers = if let Some(headers) = header_cache.get(profile_path) {
                headers.clone()
            } else {
                let headers = Self::read_matrix_headers(Path::new(profile_path))?;
                header_cache.insert(profile_path.to_string(), headers.clone());
                headers
            };

            let goldstd_headers = if let Some(headers) = header_cache.get(goldstd_path) {
                headers.clone()
            } else {
                let headers = Self::read_matrix_headers(Path::new(goldstd_path))?;
                header_cache.insert(goldstd_path.to_string(), headers.clone());
                headers
            };

            let profile_matches = match Self::match_headers_with_regex(
                &profile_headers,
                profile_regex,
                "ProfileRegex",
            ) {
                Ok(matches) => matches,
                Err(err) => {
                    if Self::match_headers_with_regex(
                        &goldstd_headers,
                        profile_regex,
                        "ProfileRegex",
                    )
                    .is_ok()
                    {
                        return Err(MetaError::DataError(format!(
                            "ProfileRegex '{}' matches GoldStd matrix columns but not Profile columns. Check for swapped regex columns.",
                            profile_regex
                        )));
                    }
                    return Err(err);
                }
            };
            let goldstd_matches = match Self::match_headers_with_regex(
                &goldstd_headers,
                goldstd_regex,
                "GoldStdRegex",
            ) {
                Ok(matches) => matches,
                Err(err) => {
                    if Self::match_headers_with_regex(
                        &profile_headers,
                        goldstd_regex,
                        "GoldStdRegex",
                    )
                    .is_ok()
                    {
                        return Err(MetaError::DataError(format!(
                            "GoldStdRegex '{}' matches Profile matrix columns but not GoldStd columns. Check for swapped regex columns.",
                            goldstd_regex
                        )));
                    }
                    return Err(err);
                }
            };

            let mut keys = profile_matches
                .keys()
                .filter(|key| goldstd_matches.contains_key(*key))
                .cloned()
                .collect::<Vec<_>>();
            keys.sort();

            if keys.is_empty() {
                return Err(MetaError::DataError(format!(
                    "No matching samples between ProfileRegex '{}' and GoldStdRegex '{}' for row {}",
                    profile_regex, goldstd_regex, row_index
                )));
            }

            for key in keys {
                let profile_sample = profile_matches.get(&key).unwrap().to_string();
                let goldstd_sample = goldstd_matches.get(&key).unwrap().to_string();

                debug!(
                    "Matrix match: Profile column '{}' <-> GoldStd column '{}' (key {})",
                    profile_sample, goldstd_sample, key
                );

                let id = format!("{}_{}", base_id, key);
                let profile_ref =
                    format!("matrix:{}::{}", profile_path, profile_sample);
                let goldstd_ref =
                    format!("matrix:{}::{}", goldstd_path, goldstd_sample);

                let entry = MetaEntry {
                    id: id.clone(),
                    sample: goldstd_sample.clone(),
                    dataset: Self::normalize_optional_str(dataset),
                    tool: tool.to_string(),
                    taxonomy: taxonomy.to_string(),
                    profile: PathBuf::from(profile_ref.clone()),
                    profile_columns: None,
                    goldstd: PathBuf::from(goldstd_ref.clone()),
                    goldstd_columns: None,
                    goldstd_tree: Self::normalize_optional_str(goldstd_tree)
                        .map(PathBuf::from),
                    taxa_list: Self::normalize_optional_str(taxa_list).map(PathBuf::from),
                };

                entries.push(entry);
                ids.push(id);
                samples.push(goldstd_sample);
                datasets.push(dataset.unwrap_or("").to_string());
                tools.push(tool.to_string());
                taxonomies.push(taxonomy.to_string());
                profiles.push(profile_ref);
                goldstds.push(goldstd_ref);
                goldstd_trees.push(goldstd_tree.unwrap_or("NA").to_string());
                taxa_lists.push(taxa_list.unwrap_or("NA").to_string());
            }
        }

        let raw = DataFrame::new(vec![
            Series::new(MetaColumnStrings::ID.into(), ids),
            Series::new(MetaColumnStrings::SAMPLE.into(), samples),
            Series::new(MetaColumnStrings::DATASET.into(), datasets),
            Series::new(MetaColumnStrings::TOOL.into(), tools),
            Series::new(MetaColumnStrings::TAXONOMY.into(), taxonomies),
            Series::new(MetaColumnStrings::PROFILE.into(), profiles),
            Series::new(MetaColumnStrings::GOLDSTD.into(), goldstds),
            Series::new(MetaColumnStrings::GOLDSTD_TREE.into(), goldstd_trees),
            Series::new(MetaColumnStrings::AVAILABLE_TAXA.into(), taxa_lists),
        ])
        .map_err(|e| MetaError::DataError(format!("Failed to build meta DataFrame: {}", e)))?;

        Ok(Self { raw, entries })
    }

    /// Loads a meta file into a Polars DataFrame based on extension.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to `.xlsx`, `.csv`, or `.tsv` meta file
    ///
    /// # Returns
    ///
    /// DataFrame if the file can be read.
    pub fn polars_from_path(path: impl AsRef<Path>) -> Option<DataFrame> {
        let ext = match path.as_ref().extension() {
            Some(ext) => ext,
            None => panic!("File needs to be either .xlsx, .csv, or .tsv. Found not extension"),
        };
        match ext.to_str() {
            Some("xlsx") => Some(Self::df_from_xlsx(path)),
            Some("tsv") => Some(Self::df_from_text(path)),
            Some("csv") => Some(Self::df_from_text(path)),
            _ => panic!("Found extension is not valid ({:?})", ext),
        }
    }

    /// Loads and parses a meta file from disk.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to `.xlsx`, `.csv`, or `.tsv` meta file
    ///
    /// # Returns
    ///
    /// Parsed `Meta` with expanded entries.
    ///
    /// # Errors
    ///
    /// Returns `MetaError` when loading or validation fails.
    pub fn from_path(path: impl AsRef<Path>) -> MetaResult {
        let ext = match path.as_ref().extension() {
            Some(ext) => ext,
            None => panic!("File needs to be either .xlsx, .csv, or .tsv. Found not extension"),
        };
        let entries = match ext.to_str() {
            Some("xlsx") => Self::df_from_xlsx(path),
            Some("tsv") => Self::df_from_text(path),
            Some("csv") => Self::df_from_text(path),
            _ => panic!("Found extension is not valid ({:?})", ext),
        };
        let is_matrix = Self::is_matrix_format(&entries);
        let meta = Self::from_polars_df(entries)?;
        if is_matrix {
            Ok(meta)
        } else {
            meta.validate()
        }
    }

    /// Returns a reference to a column by name.
    ///
    /// # Arguments
    ///
    /// * `column_name` - Column name to retrieve
    ///
    /// # Returns
    ///
    /// Column series if present.
    pub fn get_column(&self, column_name: &str) -> Option<&Series> {
        self.raw.get_columns().get(
            self.raw
                .get_column_index(column_name)
                .expect(&format!("Cannot find column '{}'", column_name)),
        )
    }

    /// Returns the unique values for a column by name.
    ///
    /// # Arguments
    ///
    /// * `column_name` - Column name to retrieve
    ///
    /// # Returns
    ///
    /// Series of unique values.
    ///
    /// # Errors
    ///
    /// Returns a Polars error when the column cannot be accessed.
    pub fn get_unique_col_values(&self, column_name: &str) -> PolarsResult<Series> {
        self.raw
            .get_columns()
            .get(
                self.raw
                    .get_column_index(column_name)
                    .expect(&format!("Cannot find column '{}'", column_name)),
            )
            .map(|x| x.unique())
            .unwrap()
    }

    /// Returns the unique values for a known meta column.
    ///
    /// # Arguments
    ///
    /// * `column` - Meta column identifier
    ///
    /// # Returns
    ///
    /// Series of unique values.
    ///
    /// # Errors
    ///
    /// Returns a Polars error when the column cannot be accessed.
    pub fn get_unique_col_values_from_col(&self, column: &MetaColumn) -> PolarsResult<Series> {
        self.raw
            .get_columns()
            .get(
                self.raw
                    .get_column_index(column.to_str())
                    .expect(&format!("Cannot find column '{}'", column.to_str())),
            )
            .map(|x| x.unique())
            .unwrap()
    }

    /// Returns the list of tools defined in the meta file.
    ///
    /// # Returns
    ///
    /// Vector of tool names, or `None` if unavailable.
    pub fn get_tools(&self) -> Option<Vec<String>> {
        Some(
            self.get_unique_col_values("Tool")
                .ok()?
                .str()
                .ok()?
                .into_iter()
                .map(|x| x.unwrap().to_string())
                .collect::<Vec<_>>(),
        )
    }

    /// Returns the list of datasets defined in the meta file.
    ///
    /// # Returns
    ///
    /// Vector of dataset names, or `None` if unavailable.
    pub fn get_datasets(&self) -> Option<Vec<String>> {
        Some(
            self.get_unique_col_values("Dataset")
                .ok()?
                .str()
                .ok()?
                .into_iter()
                .map(|x| x.unwrap().to_string())
                .collect::<Vec<_>>(),
        )
    }

    /// Returns the unique set of tree paths referenced by the meta file.
    ///
    /// # Returns
    ///
    /// Set of tree file paths.
    pub fn get_tree_path_set(&self) -> HashSet<PathBuf> {
        let col_vals = self
            .get_unique_col_values_from_col(&MetaColumn::GoldStdTree)
            .expect("Cannot get tree values from Meta polars DataFrame");

        col_vals
            .str()
            .unwrap()
            .iter()
            .map(|x| x.unwrap())
            .filter(|&s| !s.is_empty() && s != "NA")
            .map(|path_str| PathBuf::from_str(path_str).unwrap())
            .collect::<HashSet<PathBuf>>()
    }
}
