use std::{collections::{HashMap, HashSet}, fs::File, path::Path, process::exit};

use itertools::izip;

use crate::{common::{Custom, TaxonomicRank, GTDB, NCBI}, format::{Auto, Columns, CAMI}, meta::Meta, profile::{LoadProfile, Profile, ProfileWrapper}, utils::{load_file_lineages_to_hashset, load_file_to_hashset}};


#[derive(thiserror::Error, Debug, Clone)]
pub enum ProfileHandlerError {
    #[error("{0}")]
    GenericError(String),
    #[error("{0}")]
    CamiFormatError(String),
    #[error("{0}")]
    TaxonomyError(String),
}

#[derive(Default)]
pub struct ProfileHandler {
    pub meta: Meta,
    pub prediction_map: HashMap<String, ProfileWrapper>,
    pub gold_std_map: HashMap<String, ProfileWrapper>,
    pub available_taxa: HashMap<String, HashSet<String>>,
}

impl ProfileHandler {
    pub fn load_prediction_profiles(&mut self, meta: &Meta, prediction_column: &str, taxonomy_column: &str, column_format_column: &str) -> Result<(), ProfileHandlerError> {

        let profiles = meta.raw.column(prediction_column).unwrap().str().unwrap().iter();
        let taxonomy = meta.raw.column(taxonomy_column).unwrap().str().unwrap().iter();
        let column_format_series = meta.raw.column(column_format_column);

        let column_format_iter: Box<dyn Iterator<Item = Option<&str>>> = match column_format_series {
            Ok(series) => Box::new(series.str().unwrap().iter()),
            Err(_) => Box::new(std::iter::repeat(None)),
        };
        
        for (path, taxonomy, column_format) in izip!(profiles, taxonomy, column_format_iter) {
            let path = path.expect("No path");
            let profile = Self::load_profile(path, taxonomy, column_format);

            match profile {
                Some(profile) => {
                    self.prediction_map.contains_key(path)
                        .then(|| return Err::<(), ProfileHandlerError>(ProfileHandlerError::GenericError("Duplicate profile".to_owned())));
                    
                    self.prediction_map.insert(path.to_owned(), profile);
                },
                _ => (),
            }
            
        }

        Ok(())
    }

    pub fn load_gold_profiles(&mut self, meta: &Meta, gold_column: &str, taxonomy_column: &str, column_format_column: &str) -> Result<(), ProfileHandlerError> {

        let profiles = meta.raw.column(gold_column).unwrap().str().unwrap().iter();
        let taxonomy = meta.raw.column(taxonomy_column).unwrap().str().unwrap().iter();
        let column_format_series = meta.raw.column(column_format_column);

        let column_format_iter: Box<dyn Iterator<Item = Option<&str>>> = match column_format_series {
            Ok(series) => Box::new(series.str().unwrap().iter()),
            Err(_) => Box::new(std::iter::repeat(None)),
        };
        
        
        let mut data: Vec<ProfileWrapper> = Vec::default();
        
        for (path, taxonomy, column_format) in izip!(profiles, taxonomy, column_format_iter) {
            let path = path.expect("No path");
            if self.gold_std_map.contains_key(path) { continue }

            let profile = Self::load_profile(path, taxonomy, column_format);
            
            match profile {
                Some(profile) => {
                    if !self.gold_std_map.contains_key(path) {
                        self.gold_std_map.insert(path.to_owned(), profile);
                    }
                },
                _ => (),
            }
            
        }

        Ok(())
    }

    pub fn load_available_taxa(&mut self, meta: &Meta) -> std::io::Result<()> {
        for row in meta.entries.iter() {
            if let Some(taxa) = &row.taxa_list {
                if self.available_taxa.contains_key(taxa.to_str().unwrap()) { continue };
                self.available_taxa.insert(taxa.to_str().unwrap().to_owned(), load_file_lineages_to_hashset(taxa)?);
            }
        }
        Ok(())
    }

    pub fn load_profile(path: impl AsRef<Path>, taxonomy: Option<impl AsRef<str>>, column_format: Option<impl AsRef<str>>) -> Option<ProfileWrapper> {

        println!("File: {:?}", path.as_ref());
        let mut file = File::open(&path).unwrap();
        let columns = column_format.map_or(None, |str| Columns::from_format_str(str.as_ref()).ok());

        let profile = match taxonomy {
            Some(s) if s.as_ref().starts_with("GTDB") => {
                match Profile::<GTDB>::load::<Auto, _>(&mut file, columns) {
                    Ok(profile) => Some(profile.wrap()),
                    Err(e) => {
                        eprintln!("GTDB Autodetect Error: {}\nFile: {}", e, path.as_ref().display());
                        None
                    },
                }
            },
            Some(s) if s.as_ref().starts_with("NCBI") => {
                match Profile::<NCBI>::load::<CAMI, _>(&mut file, None) {
                    Ok(profile) => Some(profile.wrap()),
                    Err(e) => {
                        eprintln!("NCBI CAMI Error: {}", e);
                        None
                    },
                }
            },
            Some(_) => {
                match Profile::<Custom>::load::<Auto, _>(&mut file, columns) {
                    Ok(profile) => Some(profile.wrap()),
                    Err(e) => {
                        eprintln!("CAMI Error: {}", e);
                        None
                    },
                }
            }
            None => panic!("This should not happen"),
        };

        profile
    }

    pub fn from_meta(path: impl AsRef<Path>) -> Result<Self, ProfileHandlerError> {
        // let meta = Meta::from_path(&path).expect("Meta file not valid");
        
        let polars_df = Meta::polars_from_path(&path).expect("Meta file not valid");
        let meta = Meta::from_polars_df(polars_df).unwrap();

        // eprintln!("{:?}", newmeta.entries);

        let mut res = Self::default();
        let _ = res.load_prediction_profiles(&meta, "Profile", "Taxonomy", "ProfileColumns");
        let _ = res.load_gold_profiles(&meta, "GoldStd", "Taxonomy", "GoldStdColumns");
        let _ = res.load_available_taxa(&meta);

        for (path, profile) in  &res.prediction_map {
            let uranks = profile.unique_ranks();
            let species = profile.taxa(&TaxonomicRank::Species);
            // eprintln!("Path: {}\n\t\t{:?}", path, uranks);
            if species.is_none() {
                eprintln!("NONE {} Species: {}\n\t\t{:?} ---- {:?}", profile.taxonomy() , path, species, uranks);
            }
        }
    
        res.meta = meta;


        Ok(res)
    }
}
