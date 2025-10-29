use std::{cmp::Ordering, collections::HashMap, str::FromStr};

use itertools::Itertools;
use polars::{error::{PolarsError, PolarsResult}, frame::DataFrame, prelude::{DataType, NamedFrom}, series::Series};
use strum::{EnumIter, IntoEnumIterator};
use thiserror::Error;

use crate::{common::TaxonomicRank, profile::BC};

#[derive(Error, Debug)]
pub enum DetailedDataError {
    #[error("Column `{0}` is not available in polars df but is not optional")]
    MissingColumn(String),
    #[error("Value `{value:?}` cannot be translated into a valid item for Column `{column:?}` ")]
    WrongValue {
        column: DetailedColumn,
        value: String,
    },
    #[error("PolarsError: {0}")]
    PolarsError(#[from] PolarsError),
}

pub type DC = DetailedColumn;

#[derive(Debug, Clone, EnumIter)]
pub enum DetailedColumn {
    ID,
    Name,
    Type,
    Rank,
    PredictionAbundance,
    GoldStdAbundance,
    PredictionCount,
    GoldStdCount,
    ValidTaxon,
    DetectableTaxon,
    ClosestNeighbor,
    ClosestNeighborType,
    ClosestNeighborDistance,
    ClosestNeighborDetectable,
    ClosestNeighborAbundance,
}

impl DetailedColumn {
    pub fn to_string(&self) -> String {
        match self {
            DetailedColumn::ID => "ID".to_owned(),
            DetailedColumn::Name => "Name".to_owned(),
            DetailedColumn::Type => "Type".to_owned(),
            DetailedColumn::PredictionAbundance => "PredictionAbundance".to_owned(),
            DetailedColumn::GoldStdAbundance => "GoldStdAbundance".to_owned(),
            DetailedColumn::PredictionCount => "PredictionCount".to_owned(),
            DetailedColumn::GoldStdCount => "GoldStdCount".to_owned(),
            DetailedColumn::ValidTaxon => "ValidTaxon".to_owned(),
            DetailedColumn::DetectableTaxon => "DetectableTaxon".to_owned(),
            DetailedColumn::ClosestNeighbor => "ClosestNeighbor".to_owned(),
            DetailedColumn::ClosestNeighborType => "ClosestNeighborType".to_owned(),
            DetailedColumn::ClosestNeighborDistance => "ClosestNeighborDistance".to_owned(),
            DetailedColumn::ClosestNeighborDetectable => "ClosestNeighborDetectable".to_owned(),
            DetailedColumn::ClosestNeighborAbundance => "ClosestNeighborAbundance".to_owned(),
            DetailedColumn::Rank => "Rank".to_owned(),
        }
    }

    pub fn mandatory_columns() -> &'static [DetailedColumn] {
        &[
            Self::Name,
            Self::ID,
            Self::Type,
            Self::PredictionAbundance,
            Self::GoldStdAbundance,
            Self::Rank,
            Self::PredictionCount,
            Self::GoldStdAbundance,
            Self::GoldStdCount,
            Self::ValidTaxon,
            Self::DetectableTaxon,
        ]
    }
}

impl TryFrom<&str> for DetailedColumn {
    type Error = String;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        for ele in DetailedColumn::iter() {
            if ele.to_string() == value {
                return Ok(ele);
            }
        }

        Err(format!("Is not a valid column {}", value))
    }
}

#[derive(Clone, Debug, Default)]
pub struct DetailedDataRow {
    pub id: String,
    pub name: String,
    pub rank: TaxonomicRank,
    pub bc_type: BC,
    pub prediction_abundance: Option<f64>,
    pub goldstd_abundance: Option<f64>,
    pub prediction_count: Option<usize>,
    pub goldstd_count: Option<usize>,
    pub valid_taxon: Option<bool>,
    pub detectable_taxon: Option<bool>,
    pub closest_neigbor: Option<String>,
    pub closest_neighbor_type: Option<BC>,
    pub closest_neighbor_distance: Option<f64>,
    pub closest_neighbor_detectable: Option<bool>,
    pub closest_neighbor_abundance: Option<f64>,
}

impl DetailedDataRow {

    pub fn empty_none_or_parse_err<T: FromStr>(
        field: DetailedColumn,
        value: &str,
    ) -> Result<Option<T>, DetailedDataError> {
        if value.is_empty() {
            Ok(None)
        } else {
            Ok(Some(value.parse::<T>().map_err(|_| {
                DetailedDataError::WrongValue {
                    column: field,
                    value: value.to_owned(),
                }
            })?))
        }
    }

    pub fn set_field(
        &mut self,
        field: DetailedColumn,
        value: &str,
    ) -> Result<(), DetailedDataError> {
        match field {
            DetailedColumn::ID => self.id = value.to_owned(),
            DetailedColumn::Name => self.name = value.to_owned(),
            DetailedColumn::Type => {
                self.bc_type = BC::from_str(value).map_err(|_| DetailedDataError::WrongValue {
                    column: field,
                    value: value.to_owned(),
                })?
            }
            DetailedColumn::PredictionAbundance => {
                self.prediction_abundance = Self::empty_none_or_parse_err::<f64>(field, value)?
            }
            DetailedColumn::GoldStdAbundance => {
                self.goldstd_abundance = Self::empty_none_or_parse_err::<f64>(field, value)?
            }
            DetailedColumn::PredictionCount => {
                self.prediction_count = Self::empty_none_or_parse_err::<usize>(field, value)?
            }
            DetailedColumn::GoldStdCount => {
                self.goldstd_count = Self::empty_none_or_parse_err::<usize>(field, value)?
            }
            DetailedColumn::ValidTaxon => {
                self.valid_taxon = match value {
                    "True" => Some(true),
                    "False" => Some(false),
                    "true" => Some(true),
                    "false" => Some(false),
                    "1" => Some(true),
                    "0" => Some(false),
                    _ => {
                        return Err(DetailedDataError::WrongValue {
                            column: field,
                            value: value.to_owned(),
                        })
                    } // Assuming this returns a Result or Error.
                }
            }
            DetailedColumn::DetectableTaxon => {
                self.valid_taxon = match value.to_lowercase().as_str() {
                    "true" => Some(true),
                    "false" => Some(false),
                    "1" => Some(true),
                    "0" => Some(false),
                    "yes" => Some(true),
                    "no" => Some(false),
                    "unknown" => None,
                    _ => {
                        return Err(DetailedDataError::WrongValue {
                            column: field,
                            value: value.to_owned(),
                        })
                    } // Assuming this returns a Result or Error.
                }
            }
            DetailedColumn::ClosestNeighbor => {
                self.closest_neigbor = (!value.is_empty()).then(|| value.to_owned())
            }
            DetailedColumn::ClosestNeighborType => {
                self.closest_neighbor_type = Self::empty_none_or_parse_err::<BC>(field, value)?
            }
            DetailedColumn::ClosestNeighborDistance => {
                self.closest_neighbor_distance = Self::empty_none_or_parse_err::<f64>(field, value)?
            }
            DetailedColumn::Rank => {
                self.rank =
                    TaxonomicRank::from_str(value).map_err(|_| DetailedDataError::WrongValue {
                        column: field,
                        value: value.to_owned(),
                    })?
            }
            DetailedColumn::ClosestNeighborDetectable => {
                self.closest_neighbor_detectable = match value.to_lowercase().as_str() {
                    "true" => Some(true),
                    "false" => Some(false),
                    "1" => Some(true),
                    "0" => Some(false),
                    "yes" => Some(true),
                    "no" => Some(false),
                    _ => {
                        return Err(DetailedDataError::WrongValue {
                            column: field,
                            value: value.to_owned(),
                        })
                    } // Assuming this returns a Result or Error.
                }
            }
            DetailedColumn::ClosestNeighborAbundance => {
                self.closest_neighbor_abundance = Self::empty_none_or_parse_err::<f64>(field, value)?
            }
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Default)]
pub struct DetailedData {
    data: Vec<DetailedDataRow>,
    name_map: HashMap<String, usize>,
}

impl DetailedData {
    pub fn init_map(&mut self) {
        for (index, row) in self.data.iter().enumerate() {
            self.name_map.insert(row.name.to_string(), index);
        }
    }

    pub fn from_polars_df(df: &DataFrame) -> Result<Self, DetailedDataError> {
        let mut data = Self::default();
        data.data
            .resize_with(df.height(), || DetailedDataRow::default());

        for column in DetailedColumn::mandatory_columns() {
            let string_series = df.column(&column.to_string())?.cast(&DataType::String)?;

            let col = string_series.str().unwrap().iter().enumerate();
            for (index, value) in col {
                data.data[index].set_field(column.clone(), value.unwrap())?;
            }
        }

        data.init_map();

        Ok(data)
    }

    pub fn to_polars_df(&self) -> PolarsResult<DataFrame> {
        let bool_to_str = |b: bool| -> &str {
            if b { "true" } else { "false" }
        };

        let name_series = Series::new(DetailedColumn::Name.to_string().into(), self.data().iter().map(|r| r.name.as_str()).collect_vec());
        let id_series = Series::new(DetailedColumn::ID.to_string().into(), self.data().iter().map(|r| r.id.as_str()).collect_vec());
        let type_series = Series::new(DetailedColumn::Type.to_string().into(), self.data().iter().map(|r| r.bc_type.to_string()).collect_vec());
        let rank_series = Series::new(DetailedColumn::Rank.to_string().into(), self.data().iter().map(|r| r.rank.to_string()).collect_vec());
        let goldstd_abundance_series = Series::new(DetailedColumn::GoldStdAbundance.to_string().into(), self.data().iter().map(|r| r.goldstd_abundance.unwrap_or(0f64)).collect_vec());
        let prediction_abundance_series = Series::new(DetailedColumn::PredictionAbundance.to_string().into(), self.data().iter().map(|r| r.prediction_abundance.unwrap_or(0f64)).collect_vec());
        let goldstd_count_series = Series::new(DetailedColumn::GoldStdCount.to_string().into(), self.data().iter().map(|r| r.goldstd_count.unwrap_or(0usize) as i64).collect_vec());
        let prediction_count_series = Series::new(DetailedColumn::PredictionCount.to_string().into(), self.data().iter().map(|r| r.prediction_count.unwrap_or(0usize) as i64).collect_vec());
        let valid_taxon_series = Series::new(DetailedColumn::ValidTaxon.to_string().into(), self.data().iter().map(|r| r.valid_taxon.map_or("", bool_to_str)).collect_vec());
        let detectable_taxon_series = Series::new(DetailedColumn::DetectableTaxon.to_string().into(), self.data().iter().map(|r| r.detectable_taxon.map_or("", bool_to_str)).collect_vec());
        let closest_neighbor_series = Series::new(DetailedColumn::ClosestNeighbor.to_string().into(), self.data().iter().map(|r| r.closest_neigbor.as_ref().map(|x| x.as_str()).unwrap_or("")).collect_vec());
        let closest_neighbor_type_series = Series::new(DetailedColumn::ClosestNeighborType.to_string().into(), self.data().iter().map(|r| r.closest_neighbor_type.as_ref().map_or("".to_owned(),  |x| x.to_string())).collect_vec());
        let closest_neighbor_dist_series = Series::new(DetailedColumn::ClosestNeighborDistance.to_string().into(), self.data().iter().map(|r| r.closest_neighbor_distance.unwrap_or(0f64)).collect_vec());
        let closest_neighbor_abundance_series = Series::new(DetailedColumn::ClosestNeighborAbundance.to_string().into(), self.data().iter().map(|r| r.closest_neighbor_abundance.unwrap_or(0f64)).collect_vec());
        
        
        DataFrame::new([
            name_series,
            id_series,
            type_series,
            rank_series,
            goldstd_abundance_series,
            goldstd_count_series,
            prediction_abundance_series,
            prediction_count_series,
            valid_taxon_series,
            detectable_taxon_series,
            closest_neighbor_series,
            closest_neighbor_type_series,
            closest_neighbor_dist_series,
            closest_neighbor_abundance_series].to_vec())
    }
    

    pub fn get_mut<'a>(&'a mut self, name: &str) -> Option<&'a mut DetailedDataRow> {
        let idx: &usize = self.name_map.get(name)?;
        self.data.get_mut(*idx)
    }

    pub fn data(&self) -> &Vec<DetailedDataRow> {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut Vec<DetailedDataRow> {
        &mut self.data
    }

    pub fn sort_by_key<F, K>(&mut self, f: F)
    where
        F: FnMut(&DetailedDataRow) -> K,
        K: Ord,
    {
        self.data.sort_by_key(f);
        self.init_map();
    }

    pub fn sort_by<F>(&mut self, compare: F)
    where
    F: FnMut(&DetailedDataRow, &DetailedDataRow) -> Ordering, {
        self.data.sort_by(compare);
        self.init_map();
    }
}
