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
    path::{Path, PathBuf},
    str::FromStr,
};
use strum::{EnumIter, IntoEnumIterator};

use crate::utils::workbook_to_dataframe;

fn read_excel_to_dataframe(path: impl AsRef<Path>) -> PolarsResult<DataFrame> {
    // Open the Excel file
    let mut workbook: Xlsx<_> = calamine::open_workbook(path).expect("Cannot open Excel file");

    // Access the first sheet (assuming there's at least one)
    let range = workbook
        .worksheet_range("Sheet1")
        .expect("Cannot find sheet");

    // Prepare vectors to hold data for DataFrame columns
    let mut columns: Vec<Series> = Vec::new();
    let mut column_names: Vec<String> = Vec::new();

    // Iterate over each row and column in the sheet
    for (col_idx, row) in range.rows().enumerate() {
        let mut column_data: Vec<String> = Vec::new();

        for cell in row.iter() {
            // Extract cell content and push it to column data

            let value = match cell {
                calamine::Data::Int(i) => i.to_string(),
                calamine::Data::Float(f) => f.to_string(),
                calamine::Data::String(s) => s.to_string(),
                calamine::Data::Bool(b) => b.to_string(),
                calamine::Data::DateTime(time) => time.to_string(),
                calamine::Data::DateTimeIso(time_iso) => time_iso.to_string(),
                calamine::Data::DurationIso(duration_iso) => duration_iso.to_string(),
                calamine::Data::Error(cell_error_type) => cell_error_type.to_string(),
                calamine::Data::Empty => "".to_string(),
            };

            column_data.push(value);
        }

        // Create a Series for each column, then add to DataFrame
        let series = Series::new((&format!("column_{}", col_idx)).into(), column_data);
        columns.push(series);
        column_names.push(format!("column_{}", col_idx));
    }

    // Create the DataFrame from all columns
    DataFrame::new(columns).map_err(|e| PolarsError::ComputeError(format!("{}", e).into()))
}

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
pub struct Meta {
    pub raw: DataFrame,
    pub entries: Vec<MetaEntry>,
}

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
}

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

    pub fn from_polars_df(meta: DataFrame) -> Option<Self> {
        let mut entry = MetaEntry::default();

        let mut entries = vec![MetaEntry::default(); meta.height()];

        eprintln!("Height: {} .. {}", meta.height(), entries.len());

        meta.get_column_names()
            .iter()
            .map(|x| x.to_string())
            .for_each(|name| {
                println!("{} -> {:?}", name, MetaColumn::from_string(&name));
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

        Some(Self { raw: meta, entries })
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
        println!(
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

        let meta = Meta {
            raw: entries,
            entries: Vec::default(),
        };

        meta.validate()
    }

    pub fn get_column(&self, column_name: &str) -> Option<&Series> {
        self.raw.get_columns().get(
            self.raw
                .get_column_index(column_name)
                .expect(&format!("Cannot find column '{}'", column_name)),
        )
    }

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
