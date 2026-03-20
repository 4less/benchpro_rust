use itertools::Itertools;
use lazy_static::lazy_static;
use regex::Regex;
use std::{
    borrow::Cow,
    marker::{ConstParamTy, PhantomData},
    str::FromStr,
};
use strum::{Display, EnumCount, EnumIter, IntoEnumIterator};
use thiserror::Error;

pub const AMBIGUOUS_TAXA_REGEX_STR: &str = r"^\(([^/]+)(/[^/]+)*\)$";

lazy_static! {
    static ref AMBIGUOUS_TAXA_REGEX: Regex = Regex::new(AMBIGUOUS_TAXA_REGEX_STR).unwrap();
}

type TaxID = u64;
type TaxonMap = micromap::Map<TaxonomicRank, Taxon, { TaxonomicRank::COUNT }>;

/// Supported taxonomic ranks.
#[derive(Debug, EnumIter, EnumCount, PartialEq, PartialOrd, Eq, Ord, Clone, Hash, Display)]
pub enum TaxonomicRank {
    Superkingdom,
    Domain,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Strain,
    Unknown,
}

impl Default for TaxonomicRank {
    fn default() -> Self {
        TaxonomicRank::Unknown
    }
}

impl TaxonomicRank {
    /// Returns all supported ranks in enum order.
    ///
    /// # Returns
    ///
    /// Vector of all `TaxonomicRank` variants.
    pub fn all() -> Vec<TaxonomicRank> {
        TaxonomicRank::iter().collect_vec()
    }

    /// Parses a GTDB-style rank prefix into a `TaxonomicRank`.
    ///
    /// # Arguments
    ///
    /// * `str` - Token with GTDB prefix, e.g. `d__Bacteria`
    ///
    /// # Returns
    ///
    /// Matching rank or `None` if the prefix is not recognized.
    pub fn from_gtdb_prefix(str: &str) -> Option<Self> {
        if str.len() < 3 {
            return None;
        };
        match &str[..3] {
            "d__" => Some(Self::Domain),
            "k__" => Some(Self::Superkingdom),
            "p__" => Some(Self::Phylum),
            "c__" => Some(Self::Class),
            "o__" => Some(Self::Order),
            "f__" => Some(Self::Family),
            "g__" => Some(Self::Genus),
            "s__" => Some(Self::Species),
            "t__" => Some(Self::Strain),
            _ => None,
        }
    }

    /// Parses a rank name into a `TaxonomicRank`.
    ///
    /// # Arguments
    ///
    /// * `str` - Rank name such as `species` or `genus`
    ///
    /// # Returns
    ///
    /// Matching rank or `None` if the name is not recognized.
    pub fn from_string(str: &str) -> Option<Self> {
        match &str.to_lowercase().as_str() {
            &"domain" => Some(Self::Domain),
            &"superkingdom" => Some(Self::Superkingdom),
            &"phylum" => Some(Self::Phylum),
            &"class" => Some(Self::Class),
            &"order" => Some(Self::Order),
            &"family" => Some(Self::Family),
            &"genus" => Some(Self::Genus),
            &"species" => Some(Self::Species),
            &"strain" => Some(Self::Strain),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests_lineage {
    use super::{TaxonomicRank, Taxonomy, GTDB};

    #[test]
    fn test_lowest_rank_prefers_order_over_class() {
        let lineage = GTDB::lineage_from_string(
            "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales",
            None,
        );
        let lowest = lineage.lowest().expect("Expected lowest rank");
        assert_eq!(lowest.rank.as_ref(), Some(&TaxonomicRank::Order));
    }
}
impl FromStr for TaxonomicRank {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match &s.to_lowercase().as_str() {
            &"domain" => Self::Domain,
            &"superkingdom" => Self::Superkingdom,
            &"phylum" => Self::Phylum,
            &"class" => Self::Class,
            &"order" => Self::Order,
            &"family" => Self::Family,
            &"genus" => Self::Genus,
            &"species" => Self::Species,
            &"strain" => Self::Strain,
            _ => {
                return Err(format!(
                    "`{}` cannot be converted to a valid taxonomic rank",
                    s
                ))
            }
        })
    }
}

/// Taxonomic lineage keyed by rank.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Lineage<T: Taxonomy> {
    data: TaxonMap,
    marker: PhantomData<T>,
}

/// Parses a lineage from a string representation.
pub trait LineageFromString<T: Taxonomy> {
    /// Parses a lineage from a delimited string.
    ///
    /// # Arguments
    ///
    /// * `str` - Lineage string
    ///
    /// # Returns
    ///
    /// Parsed lineage.
    fn from_string(str: &str) -> Lineage<T>;
    /// Parses a lineage from a string using explicit rank columns.
    ///
    /// # Arguments
    ///
    /// * `str` - Lineage string
    /// * `ranks` - Ordered rank list aligned with tokens
    ///
    /// # Returns
    ///
    /// Parsed lineage.
    fn from_string_with_columns(str: &str, ranks: &Vec<TaxonomicRank>) -> Lineage<T>;
}

impl<T: Taxonomy> Lineage<T> {
    /// Returns the taxon at a specific rank.
    ///
    /// # Arguments
    ///
    /// * `rank` - Rank to retrieve
    ///
    /// # Returns
    ///
    /// Taxon at the rank or `None` if missing.
    pub fn get(&self, rank: &TaxonomicRank) -> Option<&Taxon> {
        self.data.get(rank)
    }

    /// Returns the deepest (lowest) taxon in the lineage.
    ///
    /// # Returns
    ///
    /// Lowest taxon, or `None` if empty.
    pub fn lowest(&self) -> Option<&Taxon> {
        self.data.iter().max_by_key(|x| x.0).map(|(_, taxon)| taxon)
    }

    /// Normalizes taxon names in-place for matching.
    ///
    /// # Returns
    ///
    /// `()` after mutation.
    pub fn normalize(&mut self) {
        self.data.iter_mut().for_each(|(_r, t)| t.normalize_name());
    }

    /// Returns a debug-friendly string representation.
    ///
    /// # Returns
    ///
    /// String with rank/name tuples.
    pub fn to_string(&self) -> String {
        self.data
            .iter()
            .map(|(k, v)| format!("({},{:?})", k.to_string(), v))
            .join(";")
    }
}

/// Errors produced while parsing taxonomy lineages.
#[derive(Error, Debug)]
pub enum TaxonomyLineageError {
    // #[error("data store disconnected")]
    // Disconnect(#[from] std::io::Error),f
    // #[error("the data for key `{0}` is not available")]
    // Redaction(String),
    // #[error("invalid header (expected {expected:?}, found {found:?})")]
    // InvalidHeader {
    //     expected: String,
    //     found: String,
    // },
    #[error("unknown data store error")]
    Unknown,
}

/// Enumeration of supported taxonomy systems.
#[derive(Debug, Clone, PartialEq, Eq, ConstParamTy, Display)]
pub enum TaxonomyEnum {
    NCBI,
    GTDB,
    ChocoPhlAn,
    Custom,
}

/// Trait implemented by taxonomy systems with lineage parsing.
pub trait Taxonomy {
    /// Parses a lineage string into a taxonomy-specific lineage.
    ///
    /// # Arguments
    ///
    /// * `str` - Lineage string to parse
    /// * `ranks` - Optional rank list aligned with tokens
    ///
    /// # Returns
    ///
    /// Parsed lineage.
    fn lineage_from_string(str: &str, ranks: Option<&Vec<TaxonomicRank>>) -> Lineage<Self>
    where
        Self: Sized;
    /// Returns the enum identifier for this taxonomy.
    fn get_enum() -> TaxonomyEnum;
}
/// NCBI taxonomy marker type.
#[derive(Default, Debug, Clone)]
pub struct NCBI;
/// GTDB taxonomy marker type.
#[derive(Default, Debug, Clone)]
pub struct GTDB;
/// ChocoPhlAn taxonomy marker type.
#[derive(Default, Debug, Clone)]
pub struct ChocoPhlAn;
/// Custom taxonomy marker type.
#[derive(Default, Debug, Clone)]
pub struct Custom;

pub const LINEAGE_DELIMITERS: &[char] = &[';', '|'];

impl Taxonomy for GTDB {
    // impl LineageFromString<GTDB> for Lineage<GTDB> {
    fn lineage_from_string(str: &str, _ranks: Option<&Vec<TaxonomicRank>>) -> Lineage<GTDB> {
        let tokens = str.split(|c| LINEAGE_DELIMITERS.iter().any(|&del| del == c));

        let mut res = Lineage {
            data: TaxonMap::default(),
            marker: PhantomData::<_>,
        };
        let mut last_rank: Option<TaxonomicRank> = None;

        for token in tokens {
            // println!("Token: {}", token);
            let rank = TaxonomicRank::from_gtdb_prefix(token);
            if let Some(rank) = rank {
                res.data.insert(
                    rank.clone(),
                    Taxon {
                        name: token.to_string(),
                        rank: Some(rank.clone()),
                        id: None,
                        alternative_names: None,
                    },
                );
                last_rank = Some(rank);
            } else if should_insert_species_from_unprefixed(&last_rank, &res) {
                res.data.insert(
                    TaxonomicRank::Species,
                    Taxon {
                        name: format!("s__{}", token.trim()),
                        rank: Some(TaxonomicRank::Species),
                        id: None,
                        alternative_names: None,
                    },
                );
            }
        }

        res
    }

    fn get_enum() -> TaxonomyEnum {
        return TaxonomyEnum::GTDB;
    }
}

impl Taxonomy for ChocoPhlAn {
    fn lineage_from_string(str: &str, _ranks: Option<&Vec<TaxonomicRank>>) -> Lineage<ChocoPhlAn> {
        let tokens = str.split(|c| LINEAGE_DELIMITERS.iter().any(|&del| del == c));

        let mut res = Lineage {
            data: TaxonMap::default(),
            marker: PhantomData::<_>,
        };
        let mut last_rank: Option<TaxonomicRank> = None;

        for token in tokens {
            let rank = TaxonomicRank::from_gtdb_prefix(token);
            if let Some(rank) = rank {
                res.data.insert(
                    rank.clone(),
                    Taxon {
                        name: token.to_string(),
                        rank: Some(rank.clone()),
                        id: None,
                        alternative_names: None,
                    },
                );
                last_rank = Some(rank);
            } else if should_insert_species_from_unprefixed(&last_rank, &res) {
                res.data.insert(
                    TaxonomicRank::Species,
                    Taxon {
                        name: format!("s__{}", token.trim()),
                        rank: Some(TaxonomicRank::Species),
                        id: None,
                        alternative_names: None,
                    },
                );
            }
        }

        res
    }

    fn get_enum() -> TaxonomyEnum {
        TaxonomyEnum::ChocoPhlAn
    }
}

impl Taxonomy for Custom {
    // impl LineageFromString<Custom> for  Lineage<Custom> {
    fn lineage_from_string(str: &str, _ranks: Option<&Vec<TaxonomicRank>>) -> Lineage<Self> {
        let tokens = str.split(";");

        let mut res = Lineage {
            data: TaxonMap::default(),
            marker: PhantomData::<_>,
        };
        let mut last_rank: Option<TaxonomicRank> = None;

        for token in tokens {
            let rank = TaxonomicRank::from_gtdb_prefix(token);
            if let Some(rank) = rank {
                res.data.insert(
                    rank.clone(),
                    Taxon {
                        name: token.to_string(),
                        rank: Some(rank.clone()),
                        id: None,
                        alternative_names: None,
                    },
                );
                last_rank = Some(rank);
            } else if should_insert_species_from_unprefixed(&last_rank, &res) {
                res.data.insert(
                    TaxonomicRank::Species,
                    Taxon {
                        name: format!("s__{}", token.trim()),
                        rank: Some(TaxonomicRank::Species),
                        id: None,
                        alternative_names: None,
                    },
                );
            }
        }

        res
    }

    fn get_enum() -> TaxonomyEnum {
        TaxonomyEnum::Custom
    }
}

impl Taxonomy for NCBI {
    // impl LineageFromString<NCBI> for Lineage<NCBI> {
    fn lineage_from_string(str: &str, ranks: Option<&Vec<TaxonomicRank>>) -> Lineage<Self> {
        let tokens = str.split("|");

        let mut res = Lineage {
            data: TaxonMap::default(),
            marker: PhantomData::<_>,
        };

        let mut ambiguous_tokens = None;
        for (index, token) in tokens.enumerate() {
            let token_value: Cow<'_, str> = if AMBIGUOUS_TAXA_REGEX.is_match(token) {
                let mut ambiguous_tokens_tmp = token[1..token.len() - 1]
                    .split("/")
                    .map(|s| s.to_string())
                    .collect_vec();
                let representative = ambiguous_tokens_tmp.remove(0);
                ambiguous_tokens = Some(ambiguous_tokens_tmp);
                Cow::Owned(representative)
            } else {
                Cow::Borrowed(token)
            };
            let rank = TaxonomicRank::from_gtdb_prefix(token_value.as_ref());
            if let Some(rank) = rank {
                res.data.insert(
                    rank.clone(),
                    Taxon {
                        name: token_value.to_string(),
                        rank: Some(rank),
                        id: None,
                        alternative_names: ambiguous_tokens.clone(),
                    },
                );
            }
            if let Some(ranks) = ranks.as_ref() {
                match ranks.get(index) {
                    Some(rank) => {
                        if res.data.contains_key(rank) {
                            panic!("Rank has already been defined.")
                        }
                        res.data.insert(
                            rank.clone(),
                            Taxon {
                                name: token.to_string(),
                                rank: Some(rank.clone()),
                                id: None,
                                alternative_names: ambiguous_tokens.clone(),
                            },
                        );
                    }
                    None => todo!(),
                }
            }
        }
        res.normalize();

        res
    }

    fn get_enum() -> TaxonomyEnum {
        TaxonomyEnum::NCBI
    }
}

fn should_insert_species_from_unprefixed(
    last_rank: &Option<TaxonomicRank>,
    lineage: &Lineage<impl Taxonomy>,
) -> bool {
    matches!(last_rank, Some(TaxonomicRank::Genus))
        && !lineage.data.contains_key(&TaxonomicRank::Species)
}

/// Taxon with optional rank, id, and alternative names.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct Taxon {
    pub name: String,
    pub rank: Option<TaxonomicRank>,
    pub id: Option<TaxID>,
    pub alternative_names: Option<Vec<String>>,
}

impl AsRef<Taxon> for Taxon {
    fn as_ref(&self) -> &Taxon {
        self
    }
}

impl Taxon {
    fn internal_default() -> Self {
        Self {
            name: String::default(),
            rank: None,
            id: None,
            alternative_names: None,
        }
    }

    /// Constructs a taxon with a name and rank.
    ///
    /// # Arguments
    ///
    /// * `name` - Taxon name
    /// * `rank` - Taxonomic rank
    ///
    /// # Returns
    ///
    /// Initialized `Taxon`.
    pub fn with_name_and_rank(name: &str, rank: &TaxonomicRank) -> Self {
        let mut obj = Self::internal_default();
        obj.name = name.to_string();
        obj.rank = Some(rank.to_owned());
        obj
    }

    /// Returns true if any name (primary or alternative) matches another taxon.
    ///
    /// # Arguments
    ///
    /// * `other` - Taxon to compare against
    ///
    /// # Returns
    ///
    /// True when names and ranks are compatible.
    pub fn match_any(&self, other: &Self) -> bool {
        // Early condition to throw
        if self.rank != other.rank {
            return false;
        }
        if self.name == other.name {
            return true;
        }

        return self
            .names_iter()
            .any(|self_name| other.names_iter().any(|other_name| self_name == other_name));
    }

    /// Returns true if both taxa match exactly with no alternatives.
    ///
    /// # Arguments
    ///
    /// * `other` - Taxon to compare against
    ///
    /// # Returns
    ///
    /// True when ranks and names match and no alternatives exist.
    pub fn match_exact(&self, other: &Self) -> bool {
        other.alternative_names.is_none()
            && self.alternative_names.is_none()
            && self.rank == other.rank
            && self.name == other.name
    }

    /// Iterates over the primary name and any alternative names.
    ///
    /// # Returns
    ///
    /// Iterator over taxon names.
    pub fn names_iter<'a>(&'a self) -> Box<dyn Iterator<Item = &'a str> + 'a> {
        let names = std::iter::once(self.name.as_str());
        if let Some(ambig) = &self.alternative_names {
            return Box::new(names.chain(ambig.iter().map(|name| name.as_str())));
        }
        Box::new(names)
    }

    /// Normalizes name formatting to improve matching.
    ///
    /// # Returns
    ///
    /// `()` after mutation.
    pub fn normalize_name(&mut self) {
        self.name = self.name.replace("-", "_");
        self.name = self.name.replace(" ", "_");
        self.name = self.name.replace(".", "");

        if let Some(alt) = &mut self.alternative_names {
            alt.iter_mut().for_each(|name: &mut String| {
                let new_name = name.replace("-", "_");
                let new_name = new_name.replace(" ", "_");
                name.clear();
                name.push_str(&new_name);
            });
        }
    }
}

/// Marker for whether a taxon is detectable.
pub enum Detectable {
    Unknown,
    True,
    False,
}

impl Detectable {
    /// Returns the string form used in outputs.
    ///
    /// # Returns
    ///
    /// String representation of the enum.
    pub fn to_string(&self) -> String {
        match self {
            Detectable::Unknown => "Unknown".to_string(),
            Detectable::True => "True".to_string(),
            Detectable::False => "False".to_string(),
        }
    }

    /// Parses a string into a `Detectable` enum.
    ///
    /// # Arguments
    ///
    /// * `var` - String value from outputs
    ///
    /// # Returns
    ///
    /// `Some(Detectable)` when the value is recognized.
    pub fn from_string(var: &str) -> Option<Self> {
        match var {
            "Unknown" => Some(Self::Unknown),
            "True" => Some(Self::True),
            "False" => Some(Self::False),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_lowest_rank_prefers_order_over_class() {
        let lineage = GTDB::lineage_from_string(
            "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales",
            None,
        );
        let lowest = lineage.lowest().expect("Expected lowest rank");
        assert_eq!(lowest.rank.as_ref(), Some(&TaxonomicRank::Order));
    }

    #[test]
    fn test_taxon_equality() {
        type R = TaxonomicRank;
        let a = Taxon::with_name_and_rank("Bert", &R::Species);
        let b = Taxon::with_name_and_rank("Bert", &R::Species);

        assert!(a.match_any(&b));
        assert!(b.match_any(&a));
        assert!(a.match_exact(&b));
        assert!(b.match_exact(&a));

        let a = Taxon::with_name_and_rank("Brian", &R::Species);
        let b = Taxon::with_name_and_rank("Bert", &R::Species);

        assert!(!a.match_any(&b));
        assert!(!b.match_any(&a));
        assert!(!a.match_exact(&b));
        assert!(!b.match_exact(&a));

        let a = Taxon::with_name_and_rank("Bert", &R::Class);
        let b = Taxon::with_name_and_rank("Bert", &R::Species);

        assert!(!a.match_any(&b));
        assert!(!b.match_any(&a));
        assert!(!a.match_exact(&b));
        assert!(!b.match_exact(&a));

        let mut a = Taxon::with_name_and_rank("Brian", &R::Species);
        let mut b = Taxon::with_name_and_rank("Bert", &R::Species);
        a.alternative_names = Some(vec!["Olaf".to_owned()]);
        b.alternative_names = Some(vec!["Brian".to_owned()]);

        assert!(!a.match_exact(&b));
        assert!(!b.match_exact(&a));

        assert!(a.match_any(&b));
        assert!(b.match_any(&a));

        assert!(a.match_any(&b));
        assert!(b.match_any(&a));

        a.rank = Some(R::Class);
        assert!(!a.match_any(&b));
        assert!(!b.match_any(&a));
    }

    #[test]
    fn test_load_taxonomy_gtdb() {
        type R = TaxonomicRank;
        let test = "d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter;s__Agathobacter rectalis";
        let species = "s__Agathobacter rectalis";
        let domain = "d__Bacteria";
        let phylum = "p__Bacillota_A";
        let class = "c__Clostridia";
        let order = "o__Lachnospirales";
        let family = "f__Lachnospiraceae";
        let genus = "g__Agathobacter";
        // let lineage = GTDBTaxonomy::lineage_from_string(test);
        let lineage = GTDB::lineage_from_string(test, None);

        assert_eq!(lineage.get(&R::Species).unwrap().name, species);
        assert_eq!(lineage.get(&R::Genus).unwrap().name, genus);
        assert_eq!(lineage.get(&R::Family).unwrap().name, family);
        assert_eq!(lineage.get(&R::Class).unwrap().name, class);
        assert_eq!(lineage.get(&R::Order).unwrap().name, order);
        assert_eq!(lineage.get(&R::Phylum).unwrap().name, phylum);
        assert_eq!(lineage.get(&R::Domain).unwrap().name, domain);
        assert!(matches!(lineage.get(&R::Strain), None));
    }

    #[test]
    fn test_unprefixed_after_genus_sets_species() {
        type R = TaxonomicRank;
        let test = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;\
g__Streptococcus;Incongruent [g__Streptococcus]";
        let lineage = GTDB::lineage_from_string(test, None);

        assert_eq!(
            lineage.get(&R::Species).unwrap().name,
            "s__Incongruent [g__Streptococcus]"
        );
        assert_eq!(lineage.get(&R::Genus).unwrap().name, "g__Streptococcus");
    }

    #[test]
    fn test_unprefixed_not_annotated_after_genus_sets_species() {
        type R = TaxonomicRank;
        let test = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Alteromonadales;\
f__Alteromonadaceae;g__Alteromonas;Not_annotated [g__Alteromonas]";
        let lineage = GTDB::lineage_from_string(test, None);

        assert_eq!(
            lineage.get(&R::Species).unwrap().name,
            "s__Not_annotated [g__Alteromonas]"
        );
        assert_eq!(lineage.get(&R::Genus).unwrap().name, "g__Alteromonas");
    }

    #[test]
    fn test_prefixed_species_takes_precedence() {
        type R = TaxonomicRank;
        let test = "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;\
g__Lactobacillus;s__Lactobacillus acidophilus";
        let lineage = GTDB::lineage_from_string(test, None);

        assert_eq!(
            lineage.get(&R::Species).unwrap().name,
            "s__Lactobacillus acidophilus"
        );
    }
}
