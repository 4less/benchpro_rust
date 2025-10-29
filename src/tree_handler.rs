use itertools::Itertools;
use thiserror::Error;
use std::{cell::RefCell, collections::{HashMap, HashSet}, fs, io, path::Path, process::exit};

use phylotree::tree::{NewickParseError, NodeId, Tree};

use crate::{meta::Meta, utils::{get_subtree_with_leaves, time}};




pub type TaxaSet = HashSet<String>;
pub type Name2Id = HashMap<String, NodeId>;
pub type TreeMap = HashMap<String, (Name2Id, Tree)>;

#[derive(Clone)]
pub struct TreeHandler {
    pub tree_map: TreeMap,
}

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
        tree.get_nodes().iter().chain(tree.get_leaves().iter()).for_each(|id| {
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

    pub fn from_meta(meta: &Meta) -> TreeHandlerResult<Self> {
        let mut res = TreeMap::default();

        let paths = meta.get_tree_path_set();
    

        // let res = paths.iter()
        //     .filter(|x| )
        for ele in paths.iter() {
            let key: String = ele.to_str().unwrap().to_owned();
            eprintln!("Path: {:?}", ele);
            let mut tree = Self::tree_from_file_with_cleanup(ele)?;
            Self::remove_escape_quotes(&mut tree);

            let name2id = tree.get_leaves().iter().map(|x| (*x, tree.get(x).unwrap().name.as_ref()))
                .filter(|(id, name)| name.is_some()) 
                .map(|(id, name)| (name.unwrap().to_owned(), id))
                .collect::<Name2Id>();

            // eprintln!("{:?}", name2id.keys());

            // exit(9);

            res.insert(key, (name2id, tree));
        }
        Ok(Self {
            tree_map: res
        })
    }

    pub fn get_subtree(&self, tree_path: &str, taxa: &TaxaSet) -> Option<Tree> {
        let is_valid_path = |path: &str| std::path::Path::new(path).exists();

        if !is_valid_path(tree_path) { return None };


        let (duration, (name2id, tree)) = time(|| self.tree_map.get(tree_path).unwrap());
        println!("\tGetting tree took {:?}... name2id size is : {}", duration, name2id.len());

        let (duration, ids) = time(|| taxa.iter()
            .filter(|&name| name2id.contains_key(name))
            .map(|name| { tree.get(name2id.get(name).unwrap()).unwrap().id })
            .collect_vec());

        println!("\tGetting ids took {:?} ... {}", duration, ids.len());

        
        let (duration, result) = time(|| get_subtree_with_leaves(tree, &ids, true));

        println!("\tGet Subtree with leaves took {:?}", duration);
        let result = match result {
            Ok(t) => Some(t),
            Err(e) => {
                eprintln!("Error: {}\n{:?}\n{:?}", e, ids, taxa);
                None
            },
        };
        result
    }
}

