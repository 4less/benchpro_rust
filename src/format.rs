use std::{
    fmt::Display,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    sync::atomic::{AtomicBool, Ordering},
};

use itertools::Itertools;

use crate::{
    common::{TaxonomicRank, Taxonomy},
    profile::{Entry, Profile, ProfileError, ProfileResult},
};

/// Header keywords used to identify the taxon name or id column.
pub const NAME_KEYWORDS: &[&str] = &["TAXID", "NAME", "ID"];
/// Header keywords used to identify the rank column.
pub const RANK_KEYWORDS: &[&str] = &["RANK"];
/// Header keywords used to identify the lineage column.
pub const LINEAGE_KEYWORDS: &[&str] = &["TAXPATHSN", "Lineage", "clade_name"];
/// Header keywords used to identify lineage id columns.
pub const LINEAGE_ID_KEYWORDS: &[&str] = &["TAXPATH"];
/// Header keywords used to identify abundance columns.
pub const ABUNDANCE_KEYWORDS: &[&str] = &["PERCENTAGE", "abundance", "relative_abundance"];

static CAMI_IGNORE_LINEAGE_ERROR: AtomicBool = AtomicBool::new(false);

/// Enable or disable ignoring CAMI lineage length mismatches.
///
/// # Arguments
///
/// * `value` - When true, lineage length mismatches use the last token as name.
pub fn set_cami_ignore_lineage_error(value: bool) {
    CAMI_IGNORE_LINEAGE_ERROR.store(value, Ordering::Relaxed);
}

/// Column index mapping for a taxonomic profile file.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Columns {
    pub taxon_name: Option<usize>,
    pub taxon_id: Option<usize>,
    pub rank: Option<usize>,
    pub lineage: Option<usize>,
    pub lineage_ids: Option<usize>,
    pub abundance: Option<usize>,
}

impl Display for Columns {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            r"{{\taxon_name: {:?}\n\taxon_id: {:?}\n\rank: {:?}\n\lineage: {:?}\n\lineage_ids: {:?}\n\abundance: {:?}\n}}",
            self.taxon_name,
            self.taxon_id,
            self.rank,
            self.lineage,
            self.lineage_ids,
            self.abundance
        )
    }
}

type ColumnResult = Result<Columns, ProfileError>;
impl Columns {
    /// Parses a custom format string into column indices.
    ///
    /// # Arguments
    ///
    /// * `str` - Format string in `Name|Lineage|Abundance|Rank` order
    ///
    /// # Returns
    ///
    /// Column mapping when the format is valid.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError::FormatError` for invalid format strings.
    pub fn from_format_str(str: &str) -> Result<Self, ProfileError> {
        let tokens = str.split('|').collect::<Vec<_>>();

        // column format is Name|Lineage|Abundance|Optional(Rank)
        if tokens.len() < 3 || tokens.len() > 4 {
            return Err(ProfileError::FormatError(format!("Invalid column format string ({}). Needs to be either of length 3 or 4 after splitting by '|'", str)));
        }

        let mut res = Columns::default();

        res.lineage = match tokens[1].parse::<usize>() {
            Ok(col) => Some(col),
            Err(_) => {
                return Err(ProfileError::FormatError(format!(
                    "Invalid column format string ({})",
                    str
                )))
            }
        };
        res.abundance = match tokens[2].parse::<usize>() {
            Ok(col) => Some(col),
            Err(_) => {
                return Err(ProfileError::FormatError(format!(
                    "Invalid column format string ({})",
                    str
                )))
            }
        };
        res.taxon_name = match tokens[0].parse::<usize>() {
            Ok(col) => Some(col),
            Err(_) => {
                if tokens[0] != "X" {
                    return Err(ProfileError::FormatError(format!("Invalid column format string ({}). Missing field at position 1 (Name) needs to be indicated by 'X', is {} instead", str, tokens[0])));
                }
                None
            }
        };

        if tokens.len() == 4 {
            res.rank = match tokens[3].parse::<usize>() {
                Ok(col) => Some(col),
                Err(_) => {
                    return Err(ProfileError::FormatError(format!(
                        "Invalid column format string ({})",
                        str
                    )))
                }
            };
        }

        Ok(res)
    }

    /// Returns true if a lineage string appears to be in GTDB format.
    ///
    /// # Arguments
    ///
    /// * `str` - Lineage string to inspect
    ///
    /// # Returns
    ///
    /// True when the string matches GTDB rank prefixes.
    pub fn is_gtdb_lineage(str: &str) -> bool {
        let gtdb_delimiter = ";";
        // conditions:
        // 1: splittable by ';'
        // 2: each token contains a GTDB-style rank prefix (e.g., "g__")
        let subtokens = str.split(gtdb_delimiter).collect::<Vec<_>>();

        let prefixes = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "t__"];
        subtokens
            .iter()
            .all(|token| prefixes.iter().any(|prefix| token.contains(prefix)))
    }

    /// Attempts to infer columns from a tokenized row.
    ///
    /// # Arguments
    ///
    /// * `tokens` - Tokenized row values
    ///
    /// # Returns
    ///
    /// Column mapping when inference succeeds.
    pub fn auto<T: AsRef<str>>(tokens: &[T]) -> Option<Self> {
        if tokens
            .iter()
            .all(|token| !token.as_ref().parse::<f64>().is_ok())
        {
            return None;
        }

        let abundance_column = tokens
            .iter()
            .enumerate()
            .filter(|(_index, token)| token.as_ref().parse::<f64>().is_ok())
            .collect::<Vec<_>>();

        let gtdb_lineage_column = tokens
            .iter()
            .position(|token| Self::is_gtdb_lineage(token.as_ref()));

        if let Some(lineage_index) = gtdb_lineage_column {
            if abundance_column.is_empty() {
                return None;
            }

            let numeric_indices = abundance_column.iter().map(|(index, _)| *index);
            let float_candidates = numeric_indices.clone().filter(|index| {
                let token = tokens[*index].as_ref();
                let is_numeric = token.parse::<f64>().is_ok();
                let has_decimal = token.contains('.');
                let has_exponent = token.contains('e') || token.contains('E');
                is_numeric && (has_decimal || has_exponent)
            });

            let abundance_index = float_candidates
                .filter(|index| *index != lineage_index)
                .last();

            let abundance_index = match abundance_index {
                Some(index) => index,
                None => return None,
            };

            let taxon_id_index = numeric_indices
                .filter(|index| *index != abundance_index && *index != lineage_index)
                .find(|index| {
                    let token = tokens[*index].as_ref();
                    token.parse::<i64>().is_ok()
                        && !token.contains('.')
                        && !token.contains('e')
                        && !token.contains('E')
                });

            let mut res = Columns::default();
            res.abundance = Some(abundance_index);
            res.lineage = Some(lineage_index);
            res.taxon_id = taxon_id_index;
            return Some(res);
        }

        // Fallback: if there are abundance columns, take the last one as abundance
        if !abundance_column.is_empty() {
            let mut res = Columns::default();
            res.abundance = Some(abundance_column.last().unwrap().0);
            // Find name column: first non-numeric column
            for (i, token) in tokens.iter().enumerate() {
                if token.as_ref().parse::<f64>().is_err() {
                    res.taxon_name = Some(i);
                    break;
                }
            }
            return Some(res);
        }

        None
    }

    /// Finds a CAMI header in the provided lines.
    ///
    /// # Arguments
    ///
    /// * `lines` - Candidate header lines
    ///
    /// # Returns
    ///
    /// Column mapping derived from a CAMI header.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError::CamiFormatError` if no CAMI header is found.
    pub fn find_cami_header<T: AsRef<str>>(lines: &[T]) -> ColumnResult {
        match lines
            .iter()
            .find(|&line| Self::from_cami(line.as_ref()).is_ok())
        {
            Some(str) => Self::from_cami(str.as_ref()),
            None => Err(ProfileError::CamiFormatError(
                "Unable to find CAMI header".to_owned(),
            )),
        }
    }

    /// Finds a generic header in the provided lines.
    ///
    /// # Arguments
    ///
    /// * `lines` - Candidate header lines
    ///
    /// # Returns
    ///
    /// Column mapping derived from the first valid header.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError::CamiFormatError` if no header is found.
    pub fn find_any_header<T: AsRef<str>>(lines: &[T]) -> ColumnResult {
        match lines
            .iter()
            .find(|&line| Self::from_generic_header(line.as_ref()).is_ok())
        {
            Some(str) => Self::from_generic_header(str.as_ref()),
            None => Err(ProfileError::CamiFormatError(
                "Unable to find CAMI header".to_owned(),
            )),
        }
    }

    /// Finds a MetaPhlAn header in the provided lines.
    ///
    /// # Arguments
    ///
    /// * `lines` - Candidate header lines
    ///
    /// # Returns
    ///
    /// Column mapping derived from a MetaPhlAn header.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError::FormatError` if no MetaPhlAn header is found.
    pub fn find_metaphlan_header<T: AsRef<str>>(lines: &[T]) -> ColumnResult {
        let mut has_mpa_version = false;
        let mut header_line: Option<&str> = None;

        for line in lines.iter().map(|line| line.as_ref()) {
            if line.starts_with("#mpa_") {
                has_mpa_version = true;
            }
            if line.starts_with("#clade_name") {
                header_line = Some(line);
            }
        }

        if !has_mpa_version || header_line.is_none() {
            return Err(ProfileError::FormatError(
                "Unable to find MetaPhlAn header".to_owned(),
            ));
        }

        Columns::from_generic_header(header_line.expect("Header line missing"))
    }

    /// Parses a CAMI header line into column indices.
    ///
    /// # Arguments
    ///
    /// * `str` - CAMI header line starting with `@@`
    ///
    /// # Returns
    ///
    /// Column mapping.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError::CamiFormatError` on invalid headers.
    pub fn from_cami(str: &str) -> ColumnResult {
        if !str.starts_with("@@") {
            return Err(ProfileError::CamiFormatError(
                "Cami header line expected but line does not start with '@@'".to_owned(),
            ));
        }

        let tokens: Vec<_> = str[2..].split_whitespace().collect();
        let mut columns = Columns::default();

        tokens
            .iter()
            .enumerate()
            .for_each(|(column, token)| match *token {
                "TAXID" => columns.taxon_id = Some(column),
                "RANK" => columns.rank = Some(column),
                "TAXPATH" => columns.lineage_ids = Some(column),
                "TAXPATHSN" => columns.lineage = Some(column),
                "PERCENTAGE" => columns.abundance = Some(column),
                _ => {}
            });

        if !(columns.taxon_id.is_some()
            && columns.rank.is_some()
            && columns.lineage.is_some()
            && columns.lineage_ids.is_some()
            && columns.abundance.is_some())
        {
            return Err(ProfileError::CamiFormatError(
                "Not all required columns are defined.".to_owned(),
            ));
        };

        Ok(columns)
    }

    /// Parses a generic header line into column indices.
    ///
    /// # Arguments
    ///
    /// * `str` - Header line with column names
    ///
    /// # Returns
    ///
    /// Column mapping.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError::FormatError` when header parsing fails.
    pub fn from_generic_header(str: &str) -> ColumnResult {
        // eprintln!("From generic header: \n{}", str);
        if str.len() == 0 {
            return Err(ProfileError::FormatError(
                "Header keyword matches twice".to_owned(),
            ));
        }

        let str = str.trim_start_matches('#').trim_start_matches('@');

        let tokens: Vec<_> = str.split("\t").collect();
        let mut columns = Columns::default();

        fn keymatch<T: AsRef<str>>(token: &str, keywords: &[T]) -> bool {
            keywords
                .iter()
                .map(|keyword| keyword.as_ref().to_lowercase())
                .any(|keyword| keyword == token.to_lowercase())
        }

        for (column, token) in tokens.iter().enumerate() {
            match *token {
                token if keymatch(token, NAME_KEYWORDS) => {
                    if columns.taxon_id.is_some() {
                        return Err(ProfileError::FormatError(
                            "Header keyword matches twice".to_owned(),
                        ));
                    };
                    columns.taxon_id = Some(column)
                }
                token if keymatch(token, RANK_KEYWORDS) => {
                    if columns.rank.is_some() {
                        return Err(ProfileError::FormatError(
                            "Header keyword matches twice".to_owned(),
                        ));
                    };
                    columns.rank = Some(column)
                }
                token if keymatch(token, LINEAGE_KEYWORDS) => {
                    if columns.lineage.is_some() {
                        return Err(ProfileError::FormatError(
                            "Header keyword matches twice".to_owned(),
                        ));
                    };
                    columns.lineage = Some(column)
                }
                token if keymatch(token, LINEAGE_ID_KEYWORDS) => {
                    if columns.lineage_ids.is_some() {
                        return Err(ProfileError::FormatError(
                            "Header keyword matches twice".to_owned(),
                        ));
                    };
                    columns.lineage_ids = Some(column)
                }
                token if keymatch(token, ABUNDANCE_KEYWORDS) => {
                    if columns.abundance.is_some() {
                        return Err(ProfileError::FormatError(
                            "Header keyword matches twice".to_owned(),
                        ));
                    };
                    columns.abundance = Some(column)
                }
                _ => {}
            }
        }

        if !(columns.lineage.is_some() && columns.abundance.is_some()) {
            return Err(ProfileError::FormatError(
                "Minimum required fields are Lineage and Abundance".to_owned(),
            ));
        };

        Ok(columns)
    }

    /// Returns the maximum column index referenced by this mapping.
    ///
    /// # Returns
    ///
    /// Highest column index, or `None` if no columns are set.
    pub fn max_col(&self) -> Option<usize> {
        let list = [
            self.taxon_id,
            self.taxon_name,
            self.lineage,
            self.lineage_ids,
            self.abundance,
            self.rank,
        ];
        list.iter().filter_map(|&option| option).max()
    }
}

/// Parser trait for profile formats.
pub trait Format {
    /// Loads a profile from a reader and optional column mapping.
    ///
    /// # Arguments
    ///
    /// * `input` - Reader positioned at the start of the file
    /// * `columns` - Optional column mapping override
    ///
    /// # Returns
    ///
    /// Parsed profile with taxonomy-specific entries.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError` when parsing fails.
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(
        input: &mut R,
        columns: Option<Columns>,
    ) -> ProfileResult<T>;
}

/// Provides a printable description of a format.
pub trait ProfilePrinter {
    /// Returns the user-facing format description.
    fn print_profile() -> String;
}

/// Detected profile format variants.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProfileFormat {
    CAMI,
    Custom(Columns),
    MetaPhlAn(Columns),
    Unknown,
}

// Implement specific formats
/// CAMI profile format parser.
pub struct CAMI;
/// Auto-detecting profile format parser.
pub struct Auto;
/// Minimal parser variant (reserved for future use).
pub struct Minimal;

/// Custom column mapping parser.
pub struct Custom;

impl Custom {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(
        input: &mut R,
        columns: &Columns,
    ) -> ProfileResult<T> {
        input
            .seek(SeekFrom::Start(0))
            .expect("Unable to reverse Readable input to start");

        let reader = BufReader::new(input);
        let mut result = Profile::<T>::default();

        let mut ranks: Option<Vec<TaxonomicRank>> = None;

        for line in reader.lines() {
            let line = line.expect("Cannot read line");
            // Minimal processing logic
            let tokens: Vec<_> = line.split("\t").collect();

            if line.starts_with("@Ranks:") {
                let line = line.strip_prefix("@Ranks:").unwrap().trim();
                let tokens = line.split("|").collect_vec();
                let rank_header = tokens
                    .iter()
                    .map(|token| TaxonomicRank::from_string(token))
                    .collect_vec();

                if !rank_header.iter().all(|rank| rank.is_some()) {
                    return Err(ProfileError::CamiFormatError(format!(
                        "Invalid lineage header (@Ranks:) with ({})",
                        line
                    )));
                }
                if ranks.is_some() {
                    return Err(ProfileError::CamiFormatError(format!(
                        "Lineage header (@Ranks:) defined twice ({})",
                        line
                    )));
                }

                ranks = Some(
                    rank_header
                        .into_iter()
                        .map(|rank| rank.unwrap())
                        .collect_vec(),
                );
            }

            if Auto::skip_row(&tokens, Some(&["#", "@"]), None) {
                // eprintln!("Skip row: {:?}", tokens);
                continue;
            }

            let max_col = columns.max_col().expect("No column defined");

            if tokens.len() <= max_col {
                return Err(ProfileError::FormatError(format!(
                    "Not enough columns in line. Expected {}, found {} ({:?})",
                    max_col,
                    tokens.len(),
                    tokens
                )));
            }

            let mut entry = Entry::<T>::default();

            let lineage_column = columns.lineage.expect(&format!(
                "Lineage column index undefined. (Columns {})",
                columns
            ));
            let abundance_column = columns.abundance.expect(&format!(
                "Abundance column index undefined. (Columns {})",
                columns
            ));

            let lineage = T::lineage_from_string(tokens[lineage_column], ranks.as_ref());
            entry.lineage = Some(lineage);

            entry.abundance = match tokens[abundance_column].parse::<f64>() {
                Ok(val) => val,
                Err(_) => {
                    return Err(ProfileError::CamiFormatError(format!(
                        "Expected abundance value. Field cannot be parsed to f64: {}",
                        tokens[abundance_column]
                    )))
                }
            };

            if let Some(col) = columns.taxon_id {
                entry.taxon_name = Some(tokens[col].to_owned());
            }

            if let Some(col) = columns.rank {
                entry.rank = match TaxonomicRank::from_string(tokens[col]) {
                    Some(rank) => rank,
                    None => {
                        return Err(ProfileError::CamiFormatError(format!(
                            "'{}' is not a valid taxonomic rank",
                            tokens[col]
                        )))
                    }
                };
            } else {
                let rank = entry
                    .lineage
                    .as_ref()
                    .expect("Benchpro currently requires some sort of lineage")
                    .lowest()
                    .expect("A lineage without rank information is not valid.")
                    .rank
                    .as_ref()
                    .expect("Benchpro currently requires some sort of rank information.");
                // eprintln!("Rank {:?}", entry.lineage.as_ref().unwrap().lowest().as_ref().unwrap().rank.as_ref().unwrap());
                entry.rank = rank.clone();
            }
            result.taxa.push(entry);
        }

        Ok(result)
    }
}

impl Auto {
    /// Returns true if a tokenized row should be skipped.
    ///
    /// # Arguments
    ///
    /// * `tokens` - Tokenized row values
    /// * `exclude_at_start` - Prefixes that mark a row as skippable
    /// * `exclude_if_field_contains` - Substrings that mark a row as skippable
    ///
    /// # Returns
    ///
    /// True when the row should be ignored.
    pub fn skip_row<T: AsRef<str>, U: AsRef<str>>(
        tokens: &[T],
        exclude_at_start: Option<&[U]>,
        exclude_if_field_contains: Option<&[U]>,
    ) -> bool {
        if exclude_at_start.is_some_and(|val| {
            val.iter()
                .any(|e| tokens.iter().any(|f| f.as_ref().starts_with(e.as_ref())))
        }) {
            return true;
        }
        if exclude_if_field_contains.is_some_and(|val| {
            val.iter()
                .any(|e| tokens.iter().any(|f| f.as_ref().contains(e.as_ref())))
        }) {
            return true;
        }

        // Skip row if all tokens are not parseable to f64.
        // Profiles NEED a column that contains the abundance
        tokens
            .iter()
            .all(|token| !token.as_ref().parse::<f64>().is_ok())
    }

    /// Derives column mappings by inspecting the first data row.
    ///
    /// # Arguments
    ///
    /// * `input` - Reader positioned at the start of the file
    ///
    /// # Returns
    ///
    /// Column mapping if detection succeeds.
    pub fn derive_columns<R: Read + Seek>(input: &mut R) -> Option<Columns> {
        input
            .seek(SeekFrom::Start(0))
            .expect("Failed rewinding reader to beginning of file");

        let reader = BufReader::new(input);
        let delimiter = "\t";

        // Get first line that is not a header line
        let target = reader
            .lines()
            .map(|line| {
                line.expect("Cannot read lines.")
                    .split(delimiter)
                    .map(|s| s.to_owned())
                    .collect::<Vec<_>>()
            })
            .filter(|tokens| !Self::skip_row(tokens, Some(&["#", "@", "clade_name"]), None))
            .next();

        Columns::auto(&target?)
    }

    /// Detects the profile format for the provided reader.
    ///
    /// # Arguments
    ///
    /// * `input` - Reader positioned at the start of the file
    ///
    /// # Returns
    ///
    /// Detected profile format.
    pub fn detect<R: Read + Seek>(input: &mut R) -> ProfileFormat {
        input
            .seek(SeekFrom::Start(0))
            .expect("Failed rewinding reader to beginning of file");
        let lines = BufReader::new(&mut *input)
            .lines()
            .map(|line| line.expect("Error while reading file line by line"))
            .collect::<Vec<_>>();
        input
            .seek(SeekFrom::Start(0))
            .expect("Failed rewinding reader to beginning of file");

        if Columns::find_cami_header(&lines).is_ok() {
            return ProfileFormat::CAMI;
        }
        if let Ok(columns) = Columns::find_metaphlan_header(&lines) {
            return ProfileFormat::MetaPhlAn(columns);
        }
        let columns = Self::derive_columns(input);
        match columns {
            Some(c) => ProfileFormat::Custom(c),
            None => ProfileFormat::Unknown,
        }
    }
}

impl Format for Auto {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(
        input: &mut R,
        columns: Option<Columns>,
    ) -> ProfileResult<T> {
        input
            .seek(SeekFrom::Start(0))
            .expect("Failed rewinding reader to beginning of file");

        if columns.is_some() {
            return Custom::load_profile(input, columns.as_ref().unwrap());
        }

        match Self::detect(input) {
            ProfileFormat::CAMI => CAMI::load_profile(input, None),
            ProfileFormat::Custom(columns) => Custom::load_profile(input, &columns),
            ProfileFormat::MetaPhlAn(columns) => MetaPhlAn::load_profile(input, Some(columns)),
            ProfileFormat::Unknown => {
                input
                    .seek(SeekFrom::Start(0))
                    .expect("Failed rewinding reader to beginning of file");

                let lines = BufReader::new(&mut *input)
                    .lines()
                    .map(|line| line.expect("Error while reading file line by line"))
                    .collect::<Vec<_>>();
                input
                    .seek(SeekFrom::Start(0))
                    .expect("Failed rewinding reader to beginning of file");

                let columns = Columns::find_any_header(&lines);

                if let Ok(columns) = columns {
                    return Custom::load_profile(input, &columns);
                }

                let columns = match Auto::derive_columns(input) {
                    Some(columns) => columns,
                    None => {
                        return Err(ProfileError::GenericError("Cannot auto detect".to_owned()))
                    }
                };
                let _ = input.seek(std::io::SeekFrom::Start(0));
                Custom::load_profile(input, &columns)
            }
        }
    }
}

/// MetaPhlAn profile format parser.
pub struct MetaPhlAn;

impl Format for MetaPhlAn {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(
        input: &mut R,
        columns: Option<Columns>,
    ) -> ProfileResult<T> {
        let columns = columns.ok_or_else(|| {
            ProfileError::FormatError("MetaPhlAn format requires column definitions".to_owned())
        })?;
        Custom::load_profile(input, &columns)
    }
}

impl Format for Custom {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(
        input: &mut R,
        columns: Option<Columns>,
    ) -> ProfileResult<T> {
        let columns = columns.ok_or_else(|| {
            ProfileError::FormatError("Custom format requires column definitions".to_owned())
        })?;
        Custom::load_profile(input, &columns)
    }
}

impl Format for CAMI {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(
        input: &mut R,
        _columns: Option<Columns>,
    ) -> ProfileResult<T> {
        input
            .seek(SeekFrom::Start(0))
            .expect("Failed rewinding reader to beginning of file");

        let mut result = Profile::<T>::default();
        let mut columns = None;

        let reader = BufReader::new(input);
        let mut ranks: Option<Vec<TaxonomicRank>> = None;

        for line in reader.lines() {
            let line = line.expect("Cannot read line");

            if line.is_empty() {
                continue;
            }
            if line.starts_with('#') {
                result.unstructured_meta.push(line.clone());
            } else if line.starts_with("@") {
                // Header line
                if line.starts_with("@@") {
                    columns = Some(Columns::from_cami(&line));
                } else {
                    if line.starts_with("@Ranks:") {
                        let line = line.strip_prefix("@Ranks:").unwrap().trim();
                        let tokens = line.split("|").collect_vec();
                        let rank_header = tokens
                            .iter()
                            .map(|token| TaxonomicRank::from_string(token))
                            .collect_vec();

                        if !rank_header.iter().all(|rank| rank.is_some()) {
                            return Err(ProfileError::CamiFormatError(format!(
                                "Invalid lineage header (@Ranks:) with ({})",
                                line
                            )));
                        }
                        if ranks.is_some() {
                            return Err(ProfileError::CamiFormatError(format!(
                                "Lineage header (@Ranks:) defined twice ({})",
                                line
                            )));
                        }

                        ranks = Some(
                            rank_header
                                .into_iter()
                                .map(|rank| rank.unwrap())
                                .collect_vec(),
                        );
                    }

                    // Keywords
                    let (k, v) = line.split_once(":").ok_or_else(|| {
                        ProfileError::CamiFormatError("Expected ':' in key-value line".to_owned())
                    })?;
                    result.meta.insert(k.to_owned(), v.trim().to_owned());
                }
            } else {
                let c = columns.clone().ok_or_else(|| {
                    ProfileError::CamiFormatError("Header line missing".to_owned())
                })??;

                let tokens: Vec<_> = line.split("\t").collect();
                let max_col = c.max_col().expect("No column defined");

                if tokens.len() <= max_col {
                    return Err(ProfileError::CamiFormatError(cami_row_error(
                        format!(
                            "Not enough columns in line. Expected {}, found {}",
                            max_col,
                            tokens.len()
                        ),
                        &line,
                    )));
                }

                let mut entry = Entry::<T>::default();

                let _id_column = c.taxon_id.expect("Taxon id column index undefined");
                let lineage_column = c.lineage.expect("Lineage column index undefined.");
                let abundance_column = c.abundance.expect("Abundance column index undefined");
                let rank_column = c.rank.expect("Rank column index undefined");

                let lineage = T::lineage_from_string(tokens[lineage_column], ranks.as_ref());
                entry.abundance = match tokens[abundance_column].parse::<f64>() {
                    Ok(val) => {
                        if val == 0.0 {
                            continue;
                        };
                        val
                    }
                    Err(_) => {
                        return Err(ProfileError::CamiFormatError(cami_row_error(
                            format!(
                                "Expected abundance value. Field cannot be parsed to f64: {}",
                                tokens[abundance_column]
                            ),
                            &line,
                        )))
                    }
                };

                entry.rank = match TaxonomicRank::from_string(tokens[rank_column]) {
                    Some(rank) => rank,
                    None => {
                        return Err(ProfileError::CamiFormatError(cami_row_error(
                            format!("'{}' is not a valid taxonomic rank", tokens[rank_column]),
                            &line,
                        )))
                    }
                };
                entry.lineage = Some(lineage);
                entry.taxon_name = cami_taxon_name(tokens[lineage_column], &entry.rank, &ranks)
                    .map_err(|e| ProfileError::CamiFormatError(cami_row_error(e, &line)))?;

                result.taxa.push(entry);
            }
        }

        Ok(result)
    }
}

fn cami_taxon_name(
    lineage_str: &str,
    rank: &TaxonomicRank,
    ranks: &Option<Vec<TaxonomicRank>>,
) -> Result<Option<String>, String> {
    let tokens = lineage_str.split("|").collect_vec();
    if let Some(ranks) = ranks.as_ref() {
        if let Some(rank_index) = ranks.iter().position(|r| r == rank) {
            let expected_len = rank_index + 1;
            if tokens.len() != expected_len {
                if CAMI_IGNORE_LINEAGE_ERROR.load(Ordering::Relaxed) {
                    log::warn!(
                        "CAMI lineage length mismatch ignored: expected {}, got {}; rank {:?}",
                        expected_len,
                        tokens.len(),
                        rank
                    );
                    return Ok(tokens.last().map(|token| token.to_string()));
                }
                return Err(format!(
                    "Lineage token count ({}) does not match expected length ({}) for rank {:?}",
                    tokens.len(),
                    expected_len,
                    rank
                ));
            }
            return Ok(tokens.get(rank_index).map(|token| token.to_string()));
        }
    }
    Ok(tokens.last().map(|token| token.to_string()))
}

fn cami_row_error(message: String, row: &str) -> String {
    format!("{}; row='{}'", message, row)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::sync::{Mutex, OnceLock};

    use crate::common::{TaxonomicRank, NCBI};
    use crate::format::{set_cami_ignore_lineage_error, Auto, Columns, ProfileFormat};
    use crate::profile::{LoadProfile, Profile};

    const NCBI_TEST_DIR: &str =
        concat!(env!("CARGO_MANIFEST_DIR"), "/data/test_data/profiles/gold_standard/NCBI");
    const NCBI_MOTUS_MOUSE_DIR: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/data/test_data/profiles/predictions/mOTUs3/mouse"
    );
    static CAMI_TEST_LOCK: OnceLock<Mutex<()>> = OnceLock::new();

    #[test]
    fn cami_profiles_detect_cami() {
        let _guard = cami_test_lock();
        let files =
            list_profile_files(&[Path::new(NCBI_TEST_DIR), Path::new(NCBI_MOTUS_MOUSE_DIR)]);
        assert!(
            !files.is_empty(),
            "No .txt files found in {}",
            NCBI_TEST_DIR
        );

        for path in files {
            let mut file = fs::File::open(&path).expect("Failed to open profile");
            let format = Auto::detect(&mut file);
            assert_eq!(
                format,
                ProfileFormat::CAMI,
                "Expected CAMI format for {}",
                path.display()
            );
        }
    }

    #[test]
    fn cami_lineage_length_strict_errors_when_mismatched() {
        let _guard = cami_test_lock();
        let files =
            list_profile_files(&[Path::new(NCBI_TEST_DIR), Path::new(NCBI_MOTUS_MOUSE_DIR)]);
        assert!(
            !files.is_empty(),
            "No .txt files found in {}",
            NCBI_TEST_DIR
        );

        set_cami_ignore_lineage_error(false);

        for path in files {
            let mismatch = has_lineage_length_mismatch(&path);
            let mut file = fs::File::open(&path).expect("Failed to open profile");
            let result = Profile::<NCBI>::load::<Auto, _>(&mut file, None);

            if mismatch {
                assert!(
                    result.is_err(),
                    "Expected lineage length mismatch error for {}",
                    path.display()
                );
            } else {
                assert!(
                    result.is_ok(),
                    "Unexpected error for {}: {:?}",
                    path.display(),
                    result.err()
                );
            }
        }
    }

    #[test]
    fn cami_lineage_length_ignore_allows_load() {
        let _guard = cami_test_lock();
        let files =
            list_profile_files(&[Path::new(NCBI_TEST_DIR), Path::new(NCBI_MOTUS_MOUSE_DIR)]);
        assert!(
            !files.is_empty(),
            "No .txt files found in {}",
            NCBI_TEST_DIR
        );

        with_cami_ignore(true, || {
            for path in files {
                let mut file = fs::File::open(&path).expect("Failed to open profile");
                let result = Profile::<NCBI>::load::<Auto, _>(&mut file, None);
                assert!(
                    result.is_ok(),
                    "Expected CAMI ignore lineage error to allow load for {}: {:?}",
                    path.display(),
                    result.err()
                );
            }
        });
    }

    #[test]
    fn auto_detects_gtdb_with_taxid_and_abundance() {
        let tokens = [
            "74426",
            "d__Bacteria;p__Actinomycetota;c__Coriobacteriia",
            "0.608810",
        ];

        let columns = Columns::auto(&tokens).expect("Expected GTDB auto detection");

        assert_eq!(columns.lineage, Some(1));
        assert_eq!(columns.abundance, Some(2));
        assert_eq!(columns.taxon_id, Some(0));
    }

    fn list_profile_files(dirs: &[&Path]) -> Vec<PathBuf> {
        let mut files = Vec::new();
        for dir in dirs {
            collect_files_recursive(dir, &mut files);
        }
        files
    }

    fn collect_files_recursive(dir: &Path, out: &mut Vec<PathBuf>) {
        let entries = match fs::read_dir(dir) {
            Ok(entries) => entries,
            Err(_) => return,
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                collect_files_recursive(&path, out);
            } else if path.is_file() {
                out.push(path);
            }
        }
    }

    fn has_lineage_length_mismatch(path: &Path) -> bool {
        let content = fs::read_to_string(path).expect("Failed to read profile");
        let mut ranks: Option<Vec<TaxonomicRank>> = None;
        let mut columns: Option<Columns> = None;

        for line in content.lines() {
            if line.starts_with("@Ranks:") {
                let line = line.strip_prefix("@Ranks:").unwrap().trim();
                let tokens = line.split('|').collect::<Vec<_>>();
                let mut parsed = Vec::with_capacity(tokens.len());
                for token in tokens {
                    let rank = TaxonomicRank::from_string(token).unwrap_or_else(|| {
                        panic!("Invalid rank token '{}' in {}", token, path.display())
                    });
                    parsed.push(rank);
                }
                ranks = Some(parsed);
            }
            if line.starts_with("@@") {
                columns =
                    Some(Columns::from_cami(line).unwrap_or_else(|e| {
                        panic!("Invalid CAMI header {}: {}", path.display(), e)
                    }));
            }
            if ranks.is_some() && columns.is_some() {
                break;
            }
        }

        let ranks = ranks.unwrap_or_else(|| panic!("Missing @Ranks header in {}", path.display()));
        let columns = columns.unwrap_or_else(|| panic!("Missing @@ header in {}", path.display()));

        let rank_col = columns.rank.expect("Missing rank column");
        let lineage_col = columns.lineage.expect("Missing lineage column");
        let abundance_col = columns.abundance.expect("Missing abundance column");

        for line in content.lines() {
            if line.starts_with('#') || line.starts_with('@') || line.trim().is_empty() {
                continue;
            }
            let tokens: Vec<_> = line.split('\t').collect();
            if tokens.len() <= abundance_col
                || tokens.len() <= rank_col
                || tokens.len() <= lineage_col
            {
                return true;
            }

            let abundance = tokens[abundance_col].parse::<f64>().unwrap_or(0.0);
            if abundance == 0.0 {
                continue;
            }

            let rank = match TaxonomicRank::from_string(tokens[rank_col]) {
                Some(rank) => rank,
                None => return true,
            };
            let lineage_tokens = tokens[lineage_col].split('|').collect::<Vec<_>>();
            let rank_index = match ranks.iter().position(|r| r == &rank) {
                Some(index) => index,
                None => return true,
            };
            let expected_len = rank_index + 1;
            if lineage_tokens.len() != expected_len {
                return true;
            }
        }

        false
    }

    fn cami_test_lock() -> std::sync::MutexGuard<'static, ()> {
        CAMI_TEST_LOCK
            .get_or_init(|| Mutex::new(()))
            .lock()
            .expect("Lock poisoned")
    }

    fn with_cami_ignore<F: FnOnce()>(value: bool, f: F) {
        set_cami_ignore_lineage_error(value);
        f();
        set_cami_ignore_lineage_error(false);
    }
}
