use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, Cursor},
    path::{Path, PathBuf},
    process::exit,
};

use itertools::izip;
use log::{debug, error, warn};

use crate::{
    common::{Custom, TaxonomicRank, GTDB, NCBI},
    format::{Auto, Columns, Custom as CustomFormat, CAMI},
    meta::Meta,
    profile::{LoadProfile, Profile, ProfileWrapper},
    utils::{load_file_lineages_to_hashset, load_file_to_hashset},
};

/// Errors raised while loading profiles from meta or files.
#[derive(thiserror::Error, Debug, Clone)]
pub enum ProfileHandlerError {
    #[error("{0}")]
    GenericError(String),
    #[error("{0}")]
    CamiFormatError(String),
    #[error("{0}")]
    TaxonomyError(String),
}

/// Caches loaded profiles and auxiliary data referenced by meta.
#[derive(Default)]
pub struct ProfileHandler {
    pub meta: Meta,
    pub prediction_map: HashMap<String, ProfileWrapper>,
    pub gold_std_map: HashMap<String, ProfileWrapper>,
    pub available_taxa: HashMap<String, HashSet<String>>,
}

impl ProfileHandler {
    fn split_matrix_reference(path: &Path) -> Option<(PathBuf, String)> {
        let path_str = path.to_string_lossy();
        let prefix = "matrix:";
        if !path_str.starts_with(prefix) {
            return None;
        }
        let rest = &path_str[prefix.len()..];
        let (matrix_path, column) = rest.split_once("::")?;
        Some((PathBuf::from(matrix_path), column.to_string()))
    }

    fn load_matrix_profile(
        matrix_path: &Path,
        column: &str,
        taxonomy: Option<impl AsRef<str>>,
    ) -> Option<ProfileWrapper> {
        let file = File::open(matrix_path).ok()?;
        let mut reader = BufReader::new(file);

        let mut header = String::new();
        reader.read_line(&mut header).ok()?;
        let header = header.trim_end_matches(['\n', '\r']);
        let headers = header.split('\t').collect::<Vec<_>>();
        let column_index = headers.iter().position(|&h| h == column)?;

        let mut body = String::new();
        for line in reader.lines() {
            let line = line.ok()?;
            if line.is_empty() {
                continue;
            }
            let tokens = line.split('\t').collect::<Vec<_>>();
            if tokens.len() <= column_index {
                continue;
            }
            let lineage = tokens[0];
            if lineage.is_empty() {
                continue;
            }
            let abundance_str = tokens[column_index];
            let abundance = abundance_str.parse::<f64>().unwrap_or(0.0);
            if abundance == 0.0 {
                continue;
            }
            body.push_str(lineage);
            body.push('\t');
            body.push_str(&abundance.to_string());
            body.push('\n');
        }

        let mut cursor = Cursor::new(body);
        let columns = Some(Columns {
            lineage: Some(0),
            abundance: Some(1),
            ..Columns::default()
        });

        let profile = match taxonomy {
            Some(s) if s.as_ref().starts_with("GTDB") => {
                Profile::<GTDB>::load::<Auto, _>(&mut cursor, columns).map(|profile| profile.wrap())
            }
            Some(s) if s.as_ref().starts_with("NCBI") => {
                Profile::<NCBI>::load::<Auto, _>(&mut cursor, columns).map(|profile| profile.wrap())
            }
            Some(_) => Profile::<Custom>::load::<Auto, _>(&mut cursor, columns)
                .map(|profile| profile.wrap()),
            None => return None,
        };

        match profile {
            Ok(profile) => Some(profile),
            Err(e) => {
                warn!(
                    "Matrix load error: {}\nFile: {}\nColumn: {}",
                    e,
                    matrix_path.display(),
                    column
                );
                None
            }
        }
    }
    /// Loads prediction profiles referenced in the meta table.
    ///
    /// # Arguments
    ///
    /// * `meta` - Parsed meta table
    /// * `prediction_column` - Column holding prediction profile paths
    /// * `taxonomy_column` - Column holding taxonomy identifiers
    /// * `column_format_column` - Column holding optional column format strings
    ///
    /// # Returns
    ///
    /// `Ok(())` when loading succeeds.
    ///
    /// # Errors
    ///
    /// Returns `ProfileHandlerError` on duplicate entries or parsing failures.
    pub fn load_prediction_profiles(
        &mut self,
        meta: &Meta,
        prediction_column: &str,
        taxonomy_column: &str,
        column_format_column: &str,
    ) -> Result<(), ProfileHandlerError> {
        let profiles = meta
            .raw
            .column(prediction_column)
            .unwrap()
            .str()
            .unwrap()
            .iter();
        let taxonomy = meta
            .raw
            .column(taxonomy_column)
            .unwrap()
            .str()
            .unwrap()
            .iter();
        let column_format_series = meta.raw.column(column_format_column);

        let column_format_iter: Box<dyn Iterator<Item = Option<&str>>> = match column_format_series
        {
            Ok(series) => Box::new(series.str().unwrap().iter()),
            Err(_) => Box::new(std::iter::repeat(None)),
        };

        for (path, taxonomy, column_format) in izip!(profiles, taxonomy, column_format_iter) {
            let path = path.expect("No path");
            let profile = Self::load_profile(path, taxonomy, column_format);

            match profile {
                Some(profile) => {
                    self.prediction_map.contains_key(path).then(|| {
                        return Err::<(), ProfileHandlerError>(ProfileHandlerError::GenericError(
                            "Duplicate profile".to_owned(),
                        ));
                    });

                    self.prediction_map.insert(path.to_owned(), profile);
                }
                _ => (),
            }
        }

        Ok(())
    }

    /// Loads gold standard profiles referenced in the meta table.
    ///
    /// # Arguments
    ///
    /// * `meta` - Parsed meta table
    /// * `gold_column` - Column holding gold standard profile paths
    /// * `taxonomy_column` - Column holding taxonomy identifiers
    /// * `column_format_column` - Column holding optional column format strings
    ///
    /// # Returns
    ///
    /// `Ok(())` when loading succeeds.
    ///
    /// # Errors
    ///
    /// Returns `ProfileHandlerError` on parsing failures.
    pub fn load_gold_profiles(
        &mut self,
        meta: &Meta,
        gold_column: &str,
        taxonomy_column: &str,
        column_format_column: &str,
    ) -> Result<(), ProfileHandlerError> {
        let profiles = meta.raw.column(gold_column).unwrap().str().unwrap().iter();
        let taxonomy = meta
            .raw
            .column(taxonomy_column)
            .unwrap()
            .str()
            .unwrap()
            .iter();
        let column_format_series = meta.raw.column(column_format_column);

        let column_format_iter: Box<dyn Iterator<Item = Option<&str>>> = match column_format_series
        {
            Ok(series) => Box::new(series.str().unwrap().iter()),
            Err(_) => Box::new(std::iter::repeat(None)),
        };

        let mut data: Vec<ProfileWrapper> = Vec::default();

        for (path, taxonomy, column_format) in izip!(profiles, taxonomy, column_format_iter) {
            let path = path.expect("No path");
            if self.gold_std_map.contains_key(path) {
                continue;
            }

            let profile = Self::load_profile(path, taxonomy, column_format);

            match profile {
                Some(profile) => {
                    if !self.gold_std_map.contains_key(path) {
                        self.gold_std_map.insert(path.to_owned(), profile);
                    }
                }
                _ => (),
            }
        }

        Ok(())
    }

    /// Loads available taxa lists referenced in the meta table.
    ///
    /// # Arguments
    ///
    /// * `meta` - Parsed meta table
    ///
    /// # Returns
    ///
    /// `Ok(())` when all taxa lists are loaded.
    ///
    /// # Errors
    ///
    /// Returns I/O errors when taxa files cannot be read.
    pub fn load_available_taxa(&mut self, meta: &Meta) -> std::io::Result<()> {
        for row in meta.entries.iter() {
            if let Some(taxa) = &row.taxa_list {
                if self.available_taxa.contains_key(taxa.to_str().unwrap()) {
                    continue;
                };
                self.available_taxa.insert(
                    taxa.to_str().unwrap().to_owned(),
                    load_file_lineages_to_hashset(taxa)?,
                );
            }
        }
        Ok(())
    }

    /// Loads a single profile from disk or a matrix reference.
    ///
    /// # Arguments
    ///
    /// * `path` - Path or matrix reference
    /// * `taxonomy` - Optional taxonomy identifier
    /// * `column_format` - Optional column mapping override
    ///
    /// # Returns
    ///
    /// Parsed `ProfileWrapper` if loading succeeds.
    pub fn load_profile(
        path: impl AsRef<Path>,
        taxonomy: Option<impl AsRef<str>>,
        column_format: Option<impl AsRef<str>>,
    ) -> Option<ProfileWrapper> {
        debug!("File: {:?}", path.as_ref());
        if let Some((matrix_path, column)) = Self::split_matrix_reference(path.as_ref()) {
            return Self::load_matrix_profile(&matrix_path, &column, taxonomy);
        }

        let mut file = File::open(&path).unwrap();
        let columns = column_format.map_or(None, |str| Columns::from_format_str(str.as_ref()).ok());

        let profile = if columns.is_some() {
            match taxonomy {
                Some(s) if s.as_ref().starts_with("GTDB") => {
                    match Profile::<GTDB>::load::<CustomFormat, _>(&mut file, columns) {
                        Ok(profile) => Some(profile.wrap()),
                        Err(e) => {
                            warn!(
                                "GTDB Custom Error: {}\nFile: {}",
                                e,
                                path.as_ref().display()
                            );
                            None
                        }
                    }
                }
                Some(s) if s.as_ref().starts_with("NCBI") => {
                    match Profile::<NCBI>::load::<CustomFormat, _>(&mut file, columns) {
                        Ok(profile) => Some(profile.wrap()),
                        Err(e) => {
                            warn!(
                                "NCBI Custom Error: {}\nFile: {}",
                                e,
                                path.as_ref().display()
                            );
                            None
                        }
                    }
                }
                Some(_) => match Profile::<Custom>::load::<CustomFormat, _>(&mut file, columns) {
                    Ok(profile) => Some(profile.wrap()),
                    Err(e) => {
                        warn!("Custom Error: {}\nFile: {}", e, path.as_ref().display());
                        None
                    }
                },
                None => panic!("This should not happen"),
            }
        } else {
            match taxonomy {
                Some(s) if s.as_ref().starts_with("GTDB") => {
                    match Profile::<GTDB>::load::<Auto, _>(&mut file, columns) {
                        Ok(profile) => Some(profile.wrap()),
                        Err(e) => {
                            warn!(
                                "GTDB Autodetect Error: {}\nFile: {}",
                                e,
                                path.as_ref().display()
                            );
                            None
                        }
                    }
                }
                Some(s) if s.as_ref().starts_with("NCBI") => {
                    match Profile::<NCBI>::load::<CAMI, _>(&mut file, None) {
                        Ok(profile) => Some(profile.wrap()),
                        Err(e) => {
                            warn!("NCBI CAMI Error: {}", e);
                            None
                        }
                    }
                }
                Some(_) => match Profile::<Custom>::load::<Auto, _>(&mut file, columns) {
                    Ok(profile) => Some(profile.wrap()),
                    Err(e) => {
                        warn!("CAMI Error: {}", e);
                        None
                    }
                },
                None => panic!("This should not happen"),
            }
        };

        profile
    }

    /// Builds a handler by loading profiles referenced in a meta file.
    ///
    /// # Arguments
    ///
    /// * `path` - Meta file path
    ///
    /// # Returns
    ///
    /// Initialized `ProfileHandler`.
    ///
    /// # Errors
    ///
    /// Returns `ProfileHandlerError` when meta parsing fails.
    pub fn from_meta(path: impl AsRef<Path>) -> Result<Self, ProfileHandlerError> {
        let polars_df = Meta::polars_from_path(&path).expect("Meta file not valid");
        let meta = Meta::from_polars_df(polars_df)
            .map_err(|e| ProfileHandlerError::GenericError(e.to_string()))?;

        let mut res = Self::default();
        let _ = res.load_prediction_profiles(&meta, "Profile", "Taxonomy", "ProfileColumns");
        let _ = res.load_gold_profiles(&meta, "GoldStd", "Taxonomy", "GoldStdColumns");
        let _ = res.load_available_taxa(&meta);

        for (path, profile) in &res.prediction_map {
            let uranks = profile.unique_ranks();
            let species = profile.taxa(&TaxonomicRank::Species);
            if species.is_none() {
                warn!(
                    "NONE {} Species: {}\n\t\t{:?} ---- {:?}",
                    profile.taxonomy(),
                    path,
                    species,
                    uranks
                );
            }
        }

        res.meta = meta;

        Ok(res)
    }
}
