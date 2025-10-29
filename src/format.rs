use std::{fmt::Display, io::{BufRead, BufReader, Read, Seek, SeekFrom}};

use itertools::{all, Itertools};

use crate::{common::{Taxon, TaxonomicRank, Taxonomy, NCBI}, profile::{Entry, Profile, ProfileError, ProfileResult}};


pub const NAME_KEYWORDS: &[&str] = &["TAXID", "NAME", "ID"];
pub const RANK_KEYWORDS: &[&str] = &["RANK"];
pub const LINEAGE_KEYWORDS: &[&str] = &["TAXPATHSN", "Lineage", "clade_name"];
pub const LINEAGE_ID_KEYWORDS: &[&str] = &["TAXPATH"];
pub const ABUNDANCE_KEYWORDS: &[&str] = &["PERCENTAGE", "abundance", "relative_abundance"];


#[derive(Debug, Clone, Default)]
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
        write!(f, r"{{\taxon_name: {:?}\n\taxon_id: {:?}\n\rank: {:?}\n\lineage: {:?}\n\lineage_ids: {:?}\n\abundance: {:?}\n}}", 
            self.taxon_name,
            self.taxon_id,
            self.rank,
            self.lineage,
            self.lineage_ids,
            self.abundance)
    }
}

type ColumnResult = Result<Columns, ProfileError>;
impl Columns {
    pub fn from_format_str(str: &str) -> Result<Self,ProfileError> {
        let tokens = str.split('|').collect::<Vec<_>>();

        // column format is Name|Lineage|Abundance|Optional(Rank)
        if tokens.len() < 3 || tokens.len() > 4 {
            return Err(ProfileError::FormatError(format!("Invalid column format string ({}). Needs to be either of length 3 or 4 after splitting by '|'", str)));

        }

        let mut res = Columns::default();

        res.lineage = match tokens[1].parse::<usize>() {
            Ok(col) => Some(col),
            Err(_) => return Err(ProfileError::FormatError(format!("Invalid column format string ({})", str))),
        };
        res.abundance = match tokens[2].parse::<usize>() {
            Ok(col) => Some(col),
            Err(_) => return Err(ProfileError::FormatError(format!("Invalid column format string ({})", str))),
        };
        res.taxon_name = match tokens[0].parse::<usize>() {
            Ok(col) => Some(col),
            Err(_) => {
                if tokens[0] != "X" {
                    return Err(ProfileError::FormatError(format!("Invalid column format string ({}). Missing field at position 1 (Name) needs to be indicated by 'X', is {} instead", str, tokens[0])));
                }
                None
            },
        };

        if tokens.len() == 4 {
            res.rank = match tokens[3].parse::<usize>() {
                Ok(col) => Some(col),
                Err(_) => return Err(ProfileError::FormatError(format!("Invalid column format string ({})", str))),
            };
        }
        
        Ok(res)
    }

    pub fn is_gtdb_lineage(str: &str) -> bool {
        let gtdb_delimiter = ";";
        // conditions: 
        // 1: splittable by ';'
        // 2: is __ after first char
        let subtokens = str.split(gtdb_delimiter).collect::<Vec<_>>();

        subtokens.iter().all(|&t| {
            t.len() >= 3 && &t[1..3] == "__"
        })
    }

    pub fn auto<T: AsRef<str>>(tokens: &[T]) -> Option<Self> {
        if tokens.iter().all(|token| !token.as_ref().parse::<f64>().is_ok()) {
            return None
        }

        let abundance_column = tokens.iter().enumerate()
            .filter(|(_index, token)| token.as_ref().parse::<f64>().is_ok())
            .collect::<Vec<_>>();

        let gtdb_lineage_column = tokens.iter()
            .position(|token| Self::is_gtdb_lineage(token.as_ref()));

        if gtdb_lineage_column.is_some() && abundance_column.len() == 1 {
            let mut res = Columns::default();
            res.abundance = Some(abundance_column[0].0);
            res.lineage = Some(gtdb_lineage_column.unwrap());
            return Some(res)
        }

        None
    }
    
    pub fn find_cami_header<T: AsRef<str>>(lines: &[T]) -> ColumnResult {
        match lines.iter().find(|&line| {
            Self::from_cami(line.as_ref()).is_ok()
        }) {
            Some(str) => Self::from_cami(str.as_ref()),
            None => Err(ProfileError::CamiFormatError("Unable to find CAMI header".to_owned())),
        }
    }

    pub fn find_any_header<T: AsRef<str>>(lines: &[T]) -> ColumnResult {
        match lines.iter().find(|&line| {
            Self::from_generic_header(line.as_ref()).is_ok()
        }) {
            Some(str) => Self::from_generic_header(str.as_ref()),
            None => Err(ProfileError::CamiFormatError("Unable to find CAMI header".to_owned())),
        }
    }

    pub fn from_cami(str: &str) -> ColumnResult {
        if !str.starts_with("@@") {
            return Err(ProfileError::CamiFormatError("Cami header line expected but line does not start with '@@'".to_owned()));
        }

        let tokens: Vec<_> = str[2..].split("\t").collect();
        let mut columns = Columns::default();

        tokens.iter().enumerate().for_each(|(column, token)| {
            match *token {
                "TAXID" => columns.taxon_id = Some(column),
                "RANK" => columns.rank = Some(column),
                "TAXPATH" => columns.lineage_ids = Some(column),
                "TAXPATHSN" => columns.lineage = Some(column),
                "PERCENTAGE" => columns.abundance = Some(column),
                _ => {},
            }
        });
        
        if !(columns.taxon_id.is_some() && 
                columns.rank.is_some() &&
                columns.lineage.is_some() &&
                columns.lineage_ids.is_some() &&
                columns.abundance.is_some()) {
            return Err(ProfileError::CamiFormatError("Not all required columns are defined.".to_owned()))
        };

        Ok(columns)
    }


    pub fn from_generic_header(str: &str) -> ColumnResult {
        // eprintln!("From generic header: \n{}", str);
        if str.len() == 0 {
            return Err(ProfileError::FormatError("Header keyword matches twice".to_owned()));
        }

        let str = str
            .trim_start_matches('#')
            .trim_start_matches('@');

        let tokens: Vec<_> = str.split("\t").collect();
        let mut columns = Columns::default();

        fn keymatch<T: AsRef<str>>(token: &str, keywords: &[T]) -> bool {
            keywords.iter().map(|keyword| keyword.as_ref().to_lowercase()).any(|keyword| keyword == token.to_lowercase())
        }

        for (column, token) in tokens.iter().enumerate() {
            match *token {
                token if keymatch(token, NAME_KEYWORDS) => {
                    if columns.taxon_id.is_some() { 
                        return Err(ProfileError::FormatError("Header keyword matches twice".to_owned()));
                    };
                    columns.taxon_id = Some(column)
                },
                token if keymatch(token, RANK_KEYWORDS) => {
                    if columns.rank.is_some() { 
                        return Err(ProfileError::FormatError("Header keyword matches twice".to_owned()));
                    };
                    columns.rank = Some(column)
                },
                token if keymatch(token, LINEAGE_KEYWORDS) => {
                    if columns.lineage.is_some() { 
                        return Err(ProfileError::FormatError("Header keyword matches twice".to_owned()));
                    };
                    columns.lineage = Some(column)
                },
                token if keymatch(token, LINEAGE_ID_KEYWORDS) => {
                    if columns.lineage_ids.is_some() { 
                        return Err(ProfileError::FormatError("Header keyword matches twice".to_owned()));
                    };
                    columns.lineage_ids = Some(column)
                },
                token if keymatch(token, ABUNDANCE_KEYWORDS) => {
                    if columns.abundance.is_some() { 
                        return Err(ProfileError::FormatError("Header keyword matches twice".to_owned()));
                    };
                    columns.abundance = Some(column)
                },
                _ => {},
            }
        }
        
        if !(columns.lineage.is_some() && 
                columns.abundance.is_some()) {
            return Err(ProfileError::FormatError("Minimum required fields are Lineage and Abundance".to_owned()))
        };

        Ok(columns)
    }

    pub fn max_col(&self) -> Option<usize> {
        let list = [self.taxon_id, self.taxon_name, self.lineage, self.lineage_ids, self.abundance, self.rank];
        list.iter()
            .filter_map(|&option| option)
            .max()
    }
}


// Define the Format trait
pub trait Format {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(input: &mut R, columns: Option<Columns>) -> ProfileResult<T>;
}

pub trait ProfilePrinter {
    fn print_profile() -> String;
}



pub enum ProfileFormat {
    CAMI,
    Custom(Columns),
    Unknown,
}


// Implement specific formats
pub struct CAMI;
pub struct Auto;
pub struct Minimal;

pub struct Custom;

impl Custom {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(input: &mut R, columns: &Columns) -> ProfileResult<T> {
        input.seek(SeekFrom::Start(0)).expect("Unable to reverse Readable input to start");

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
                let rank_header = tokens.iter().map(|token| TaxonomicRank::from_string(token)).collect_vec();
                
                if !rank_header.iter().all(|rank| rank.is_some()) {
                    return Err(ProfileError::CamiFormatError(format!("Invalid lineage header (@Ranks:) with ({})", line)))
                }
                if ranks.is_some() {
                    return Err(ProfileError::CamiFormatError(format!("Lineage header (@Ranks:) defined twice ({})", line)))
                }

                ranks = Some(rank_header.into_iter().map(|rank| rank.unwrap()).collect_vec());
            }
            
            if Auto::skip_row(&tokens, Some(&["#", "@"]), None) {
                // eprintln!("Skip row: {:?}", tokens);
                continue
            }

            let max_col = columns.max_col().expect("No column defined");

            if tokens.len() <= max_col {
                return Err(ProfileError::FormatError(format!(
                    "Not enough columns in line. Expected {}, found {} ({:?})",
                    max_col, tokens.len(), tokens
                )));
            }

            let mut entry = Entry::<T>::default();

            let lineage_column = columns.lineage.expect(&format!("Lineage column index undefined. (Columns {})", columns));
            let abundance_column = columns.abundance.expect(&format!("Abundance column index undefined. (Columns {})", columns));

            let lineage = T::lineage_from_string(tokens[lineage_column], ranks.as_ref());
            entry.lineage = Some(lineage);

            entry.abundance = match tokens[abundance_column].parse::<f64>() {
                Ok(val) => val,
                Err(_) => return Err(ProfileError::CamiFormatError(format!("Expected abundance value. Field cannot be parsed to f64: {}", tokens[abundance_column]))),
            };

            if let Some(col) = columns.taxon_id {
                entry.taxon_name = Some(tokens[col].to_owned());
            }

            if let Some(col) = columns.rank {
                entry.rank = match TaxonomicRank::from_string(tokens[col]) {
                    Some(rank) => rank,
                    None => return Err(ProfileError::CamiFormatError(format!("'{}' is not a valid taxonomic rank", tokens[col]))),
                };
            } else {
                let rank = entry.lineage.as_ref().expect("Benchpro currently requires some sort of lineage").lowest().expect("A lineage without rank information is not valid.").rank.as_ref().expect("Benchpro currently requires some sort of rank information.");
                // eprintln!("Rank {:?}", entry.lineage.as_ref().unwrap().lowest().as_ref().unwrap().rank.as_ref().unwrap());
                entry.rank = rank.clone();
            }
            result.taxa.push(entry);
        }

        Ok(result)
    }
}


impl Auto {
    pub fn skip_row<T: AsRef<str>, U: AsRef<str>>(tokens: &[T], exclude_at_start: Option<&[U]>, exclude_if_field_contains: Option<&[U]>) -> bool {
        if exclude_at_start.is_some_and(|val| {
            val.iter().any(|e| {
                tokens.iter().any(|f| f.as_ref().starts_with(e.as_ref()))
            })
        }) {
            return true
        }
        if exclude_if_field_contains.is_some_and(|val| {
            val.iter().any(|e| {
                tokens.iter().any(|f| f.as_ref().contains(e.as_ref()))
            })
        }) {
            return true
        }

        // Skip row if all tokens are not parseable to f64. 
        // Profiles NEED a column that contains the abundance
        tokens.iter().all(|token| !token.as_ref().parse::<f64>().is_ok())
    }

    pub fn derive_columns<R: Read + Seek>(input: &mut R) -> Option<Columns> {
        input.seek(SeekFrom::Start(0)).expect("Failed rewinding reader to beginning of file");

        let reader = BufReader::new(input);
        let delimiter = "\t";

        // Get first line that is not a header line
        let target = reader.lines().map(|line| {
            line.expect("Cannot read lines.").split(delimiter).map(|s| s.to_owned()).collect::<Vec<_>>()
        }).filter(|tokens| !Self::skip_row(tokens, Some(&["#", "@", "clade_name"]), None))
            .next();

        Columns::auto(&target?)
    }

    pub fn detect<R: Read + Seek>(input: &mut R) -> ProfileFormat {
        input.seek(SeekFrom::Start(0)).expect("Failed rewinding reader to beginning of file");
        // Check CAMI

        let cami = CAMI::load_profile::<NCBI, _>(input, None);
        if cami.is_ok() {
            return ProfileFormat::CAMI
        }
        let columns = Self::derive_columns(input);
        match columns {
            Some(c) => ProfileFormat::Custom(c),
            None => ProfileFormat::Unknown,
        }
    }
}

impl Format for Auto {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(input: &mut R, columns: Option<Columns>) -> ProfileResult<T> {
        input.seek(SeekFrom::Start(0)).expect("Failed rewinding reader to beginning of file");

        if columns.is_some() {
            return Custom::load_profile(input, columns.as_ref().unwrap());
        }

        match Self::detect(input) {
            ProfileFormat::CAMI => CAMI::load_profile(input, None),
            ProfileFormat::Custom(columns) => Custom::load_profile(input, &columns),
            ProfileFormat::Unknown => {
                input.seek(SeekFrom::Start(0)).expect("Failed rewinding reader to beginning of file");

                let lines = BufReader::new(&mut *input)
                    .lines()
                    .map(|line| line.expect("Error while reading file line by line"))
                    .collect::<Vec<_>>();
                input.seek(SeekFrom::Start(0)).expect("Failed rewinding reader to beginning of file");

                let columns = Columns::find_any_header(&lines);
                
                if let Ok(columns) = columns {
                    return Custom::load_profile(input, &columns)
                }

                let columns = match Auto::derive_columns(input) {
                    Some(columns) => columns,
                    None => return Err(ProfileError::GenericError("Cannot auto detect".to_owned())),
                };
                let _ = input.seek(std::io::SeekFrom::Start(0));
                Custom::load_profile(input, &columns)
            },
        }
    }
}

impl Format for CAMI {
    fn load_profile<T: Taxonomy + Default, R: Read + Seek>(input: &mut R, _columns: Option<Columns>) -> ProfileResult<T> {
        input.seek(SeekFrom::Start(0)).expect("Failed rewinding reader to beginning of file");

        let mut result = Profile::<T>::default();
        let mut columns = None;

        let reader = BufReader::new(input);
        let mut ranks: Option<Vec<TaxonomicRank>> = None;

        for line in reader.lines() {
            let line = line.expect("Cannot read line");

            if line.is_empty() { continue }
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
                        let rank_header = tokens.iter().map(|token| TaxonomicRank::from_string(token)).collect_vec();
                        
                        if !rank_header.iter().all(|rank| rank.is_some()) {
                            return Err(ProfileError::CamiFormatError(format!("Invalid lineage header (@Ranks:) with ({})", line)))
                        }
                        if ranks.is_some() {
                            return Err(ProfileError::CamiFormatError(format!("Lineage header (@Ranks:) defined twice ({})", line)))
                        }
        
                        ranks = Some(rank_header.into_iter().map(|rank| rank.unwrap()).collect_vec());
                    }

                    // Keywords
                    let (k,v) = line.split_once(":")
                        .ok_or_else(|| ProfileError::CamiFormatError("Expected ':' in key-value line".to_owned()))?;
                    result.meta.insert(k.to_owned(), v.trim().to_owned());
                }
            } else {
                let c = columns.clone()
                    .ok_or_else(|| ProfileError::CamiFormatError("Header line missing".to_owned()))??;

                let tokens: Vec<_> = line.split("\t").collect();
                let max_col = c.max_col().expect("No column defined");

                if tokens.len() <= max_col {
                    return Err(ProfileError::CamiFormatError(format!(
                        "Not enough columns in line. Expected {}, found {}",
                        max_col, tokens.len()
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
                        if val == 0.0{ continue };
                        val
                    },
                    Err(_) => return Err(ProfileError::CamiFormatError(format!("Expected abundance value. Field cannot be parsed to f64: {}", tokens[abundance_column]))),
                };
                
                entry.rank = match TaxonomicRank::from_string(tokens[rank_column]) {
                    Some(rank) => rank,
                    None => return Err(ProfileError::CamiFormatError(format!("'{}' is not a valid taxonomic rank", tokens[rank_column]))),
                };
                entry.lineage = Some(lineage);

                result.taxa.push(entry);
            }
        }

        Ok(result)
    }
}
