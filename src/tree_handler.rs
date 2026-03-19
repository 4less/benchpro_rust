use itertools::Itertools;
use std::{
    collections::{HashMap, HashSet},
    fs, io,
    path::Path,
};

use phylotree::tree::{NewickParseError, NodeId, Tree};

use crate::{
    meta::Meta,
    utils::{get_subtree_with_leaves, time},
};
use log::{debug, error};

pub type TaxaSet = HashSet<String>;
pub type Name2Id = HashMap<String, NodeId>;
pub type TreeMap = HashMap<String, (Name2Id, Tree)>;

/// Caches trees and name-to-id maps for fast subtree queries.
#[derive(Clone)]
pub struct TreeHandler {
    pub tree_map: TreeMap,
}

/// Errors encountered while loading or parsing trees.
#[derive(thiserror::Error, Debug)]
pub enum TreeHandlerError {
    #[error("Newick Parse Error: {0}")]
    NewickParseError(#[from] NewickParseError),
    #[error("IO Error: {0}")]
    IOError(#[from] io::Error),
}

type TreeHandlerResult<T> = Result<T, TreeHandlerError>;

impl TreeHandler {
    fn clean_newick_str(newick_str: &str) -> String {
        let single_quotes = newick_str.chars().filter(|c| *c == '\'').count();
        let double_quotes = newick_str.chars().filter(|c| *c == '"').count();

        let mut result_newick_str = newick_str.to_owned();
        if single_quotes > 0 && double_quotes == 0 {
            result_newick_str = newick_str.replace("'", "\"");
        }

        result_newick_str
    }

    fn remove_escape_quotes(tree: &mut Tree) {
        tree.search_nodes(|_| true)
            .iter()
            .chain(tree.get_leaves().iter())
            .for_each(|id| {
                let node = tree.get_mut(id).unwrap();
                if let Some(name) = &mut node.name {
                    *name = name.replace("\"", "");
                }
            });
    }

    fn tree_from_file_with_cleanup(path: impl AsRef<Path>) -> TreeHandlerResult<Tree> {
        let raw_newick: String = fs::read_to_string(path.as_ref())?;

        let clean_newick = Self::clean_newick_str(&raw_newick);

        Ok(Tree::from_newick(&clean_newick)?)
    }

    /// Loads all tree files referenced by the meta table.
    ///
    /// # Arguments
    ///
    /// * `meta` - Parsed meta table
    ///
    /// # Returns
    ///
    /// Initialized `TreeHandler` with cached trees.
    ///
    /// # Errors
    ///
    /// Returns `TreeHandlerError` when trees cannot be read or parsed.
    pub fn from_meta(meta: &Meta) -> TreeHandlerResult<Self> {
        let mut res = TreeMap::default();

        let paths = meta.get_tree_path_set();

        for ele in paths.iter() {
            let key: String = ele.to_str().unwrap().to_owned();
            debug!("Path: {:?}", ele);
            let mut tree = Self::tree_from_file_with_cleanup(ele)?;
            Self::remove_escape_quotes(&mut tree);

            let name2id = tree
                .get_leaves()
                .iter()
                .map(|x| (*x, tree.get(x).unwrap().name.as_ref()))
                .filter(|(_id, name)| name.is_some())
                .map(|(id, name)| (name.unwrap().to_owned(), id))
                .collect::<Name2Id>();

            res.insert(key, (name2id, tree));
        }
        Ok(Self { tree_map: res })
    }

    /// Builds a subtree containing the requested taxa from a cached tree.
    ///
    /// # Arguments
    ///
    /// * `tree_path` - Path identifying the cached tree
    /// * `taxa` - Set of taxa names to include
    ///
    /// # Returns
    ///
    /// Subtree if the tree is found and taxa can be resolved.
    pub fn get_subtree(&self, tree_path: &str, taxa: &TaxaSet) -> Option<Tree> {
        let is_valid_path = |path: &str| std::path::Path::new(path).exists();

        if !is_valid_path(tree_path) {
            return None;
        };

        let (duration, (name2id, tree)) = time(|| self.tree_map.get(tree_path).unwrap());
        debug!(
            "\tGetting tree took {:?}... name2id size is : {}",
            duration,
            name2id.len()
        );

        let (duration, ids) = time(|| {
            taxa.iter()
                .filter(|&name| name2id.contains_key(name))
                .map(|name| tree.get(name2id.get(name).unwrap()).unwrap().id)
                .collect_vec()
        });

        debug!("\tGetting ids took {:?} ... {}", duration, ids.len());

        let (duration, result) = time(|| get_subtree_with_leaves(tree, &ids, true));

        debug!("\tGet Subtree with leaves took {:?}", duration);
        let result = match result {
            Ok(t) => Some(t),
            Err(e) => {
                error!("Error: {}\n{:?}\n{:?}", e, ids, taxa);
                None
            }
        };
        result
    }
}
