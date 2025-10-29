use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, Read, Seek},
    mem::swap,
    path::Path,
    time::{Duration, Instant},
};

use calamine::{Data, DataType, Reader, Xlsx};
use itertools::Itertools;
use phylotree::tree::{Edge, Node, NodeError, NodeId, Tree, TreeError};
use polars::prelude::*;

use crate::{common::LINEAGE_DELIMITERS, meta::MetaColumn};

pub fn workbook_to_dataframe(workbook: &mut Xlsx<impl Read + Seek>) -> PolarsResult<DataFrame> {
    // Access the first sheet in the workbook
    let binding = workbook.sheet_names();

    let sheet_name = binding
        .first()
        .ok_or_else(|| PolarsError::ComputeError("No sheets found in workbook".into()))?;

    let range = workbook.worksheet_range(sheet_name).expect("Cannot parse");

    // let mut columns: Vec<Series> = Vec::new();
    // let mut columns: Vec<Vec<String>> = Vec::new();
    let mut columns: Vec<Vec<Data>> = Vec::new();
    let mut column_names: Vec<String> = Vec::new();
    let mut initialized = false;

    for (row_index, row) in range.rows().enumerate() {
        if row_index == 0 {
            // First row is treated as headers
            column_names = row
                .iter()
                .map(|cell| match cell {
                    calamine::Data::String(s) => s.clone(),
                    _ => "Unknown".to_string(),
                })
                .collect();
            continue;
        }

        if !initialized {
            columns = vec![Vec::default(); row.len()];
            initialized = true;
        }

        // Collect each cell into the appropriate column vector
        for (i, cell) in row.iter().enumerate() {
            columns[i].push(cell.clone())
        }
    }

    let mut df = DataFrame::default();
    for (i, column) in columns.iter().enumerate() {
        match column.first().as_ref().unwrap() {
            Data::Int(_) => {
                let v = column
                    .iter()
                    .map(|x| x.as_i64().expect(&format!("Cannot convert to i64 ({})", x)))
                    .collect::<Vec<i64>>();

                let s = Series::new(column_names[i].clone().into(), v);
                df.with_column(s).expect("Cannot append column");
            }
            Data::Float(_) => {
                let v = column
                    .iter()
                    .map(|x| x.as_f64().expect(&format!("Cannot convert to f64 ({})", x)))
                    .collect::<Vec<f64>>();

                let s = Series::new(column_names[i].clone().into(), v);
                df.with_column(s).expect("Cannot append column");
            }
            Data::String(_) => {
                let v = column
                    .iter()
                    .map(|x| x.as_string().unwrap_or(x.to_string()))
                    .collect::<Vec<String>>();

                let s = Series::new(column_names[i].clone().into(), v);
                df.with_column(s).expect("Cannot append column");
            }
            _ => panic!("Not yet implemented"),
        }
    }
    Ok(df)
}

pub fn time<F, T>(f: F) -> (Duration, T)
where
    F: FnOnce() -> T,
{
    let start: Instant = Instant::now();
    let result = f();
    (start.elapsed(), result)
}

// Macro to flatten nested tuples
#[macro_export]
macro_rules! flatten_tuple {
    // Base case: A single tuple of elements
    ($a:expr, $b:expr) => {
        ($a, $b)
    };

    // Recursive case: Nested tuple
    (($($t:tt),*), $($rest:tt),*) => {
        let mut flattened = flatten_tuple!($($t),*);
        flattened.push(flatten_tuple!($($rest),*));
        flattened
    };
}

pub fn f1_score(tp: usize, fp: usize, fn_: usize) -> f64 {
    (tp * 2) as f64 / (2 * tp + fp + fn_) as f64
}

pub fn sensitivity(tp: usize, fp: usize, fn_: usize) -> f64 {
    tp as f64 / (tp + fn_) as f64
}

pub fn precision(tp: usize, fp: usize, fn_: usize) -> f64 {
    tp as f64 / (tp + fp) as f64
}

/// Reads a file line by line and loads each unique line into a HashSet of Strings.
/// Returns a Result with the HashSet or an error if something goes wrong.
pub fn load_file_to_hashset<P>(file_path: P) -> std::io::Result<HashSet<String>>
where
    P: AsRef<Path>,
{
    // Open the file
    let file = File::open(file_path)?;
    let reader = std::io::BufReader::new(file);

    // Initialize an empty HashSet
    let mut lines_set = HashSet::new();

    // Read the file line by line and insert into the HashSet
    for line in reader.lines() {
        let line = line?; // Handle potential I/O errors
        lines_set.insert(line);
    }

    Ok(lines_set)
}

/// Reads a file line by line and loads each unique line into a HashSet of Strings.
/// Returns a Result with the HashSet or an error if something goes wrong.
pub fn load_file_lineages_to_hashset<P>(file_path: P) -> std::io::Result<HashSet<String>>
where
    P: AsRef<Path>,
{
    // Open the file
    let file = File::open(file_path)?;
    let reader = std::io::BufReader::new(file);

    // Initialize an empty HashSet
    let mut lines_set = HashSet::new();

    // Read the file line by line and insert into the HashSet
    for line in reader.lines() {
        let line = line?; // Handle potential I/O errors

        let tokens = line
            .split(|c| LINEAGE_DELIMITERS.iter().any(|&del| del == c))
            .collect_vec();

        if tokens.len() > 1 {
            tokens.iter().for_each(|token| {
                lines_set.insert(token.to_string());
            });
        } else {
            lines_set.insert(line);
        }
    }

    Ok(lines_set)
}

pub fn write_df(
    df: &mut DataFrame,
    path: impl AsRef<Path>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Open the file for writing
    let mut file = File::create(path).expect("Could not create file");

    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(df)?;

    Ok(())
}

pub type NameValuePair = (String, String);
pub fn add_string_columns(df: &mut DataFrame, columns: &[NameValuePair]) -> PolarsResult<()> {
    for (name, value) in columns {
        df.with_column(Series::new(
            name.into(),
            std::iter::repeat_n(value, df.height())
                .map(|x| x.clone())
                .collect::<Vec<_>>(),
        ))?;
    }
    Ok(())
}

// pub fn add_columns<T, V, U>(df: &mut DataFrame, name: &str, value: T) -> PolarsResult<()>
// where
//     T: Clone + IntoSeries + 'static, // T must implement Clone and IntoSeries
//     V: NamedFrom<&'static [T], U>,          // V must be Vec<T> or something that can be treated as a slice of T
// {
//     // Create a Vec<T> by repeating the value
//     let values: V = std::iter::repeat(value.clone()).take(df.height()).collect::<Vec<T>>().into();

//     // Use AsRef to convert values into a slice and then create a Series from it
//     let series = Series::new(name.into(), values);

//     // Add the column to the DataFrame
//     df.with_column(series)?;

//     Ok(())
// }


// pub fn add_columns<T>(df: &mut DataFrame, name: &str, value: T) -> PolarsResult<()>
// where
//     T: Clone + IntoSeries + 'static, // T must implement Clone and IntoSeries
// {
//     // Create a Vec<T> by repeating the value
//     let values: Vec<T> = std::iter::repeat(value.clone()).take(df.height()).collect();

//     // Convert the Vec<T> into a Series
//     let series = Series::new(name.into(), values);

//     // Add the column to the DataFrame
//     df.with_column(series)?;

//     Ok(())
// }

pub fn sample_apply<F>(df: &DataFrame, mut f: F) -> PolarsResult<DataFrame>
where
    F: FnMut(DataFrame) -> PolarsResult<DataFrame> + Send + Sync,
{
    let result = df
        .clone()
        .group_by(["ID", "Rank", "AllowAlternatives"])
        .unwrap()
        .apply(f);

    result
}

pub fn get_node_ids_for_leaf(tree: &Tree, leaf_id: &usize) -> HashSet<usize> {
    let mut node = tree.get(leaf_id).unwrap();

    let mut res = HashSet::new();
    res.insert(*leaf_id);

    while !node.is_root() {
        let parent_id = &node
            .parent
            .expect("Each node needs a parent except root node.");
        let parent = tree.get(parent_id).unwrap();
        res.insert(*parent_id);
        node = parent;
    }
    res
}

pub fn get_node_ids_for_leaf_set(
    tree: &Tree,
    leaf_ids: &[usize],
    from_lca: bool,
) -> HashSet<usize> {
    let mut res = HashSet::new();

    let end_node = if from_lca {
        get_lca(tree, leaf_ids).unwrap()
    } else {
        tree.get_root().unwrap()
    };

    for &leaf_id in leaf_ids {
        let mut node = tree.get(&leaf_id).unwrap();

        res.insert(leaf_id);

        while node.id != end_node {
            let parent_id = &node
                .parent
                .expect("Each node needs a parent except root node.");
            let parent = tree.get(parent_id).unwrap();
            res.insert(*parent_id);
            node = parent;
        }
    }
    res
}

pub fn get_subtree(tree: &Tree, id: &usize) -> Result<Tree, TreeError> {
    let mut new_tree = Tree::new();

    let mut new_root = tree.get(id)?.clone();

    new_root.parent = None;

    add_recursively(tree, &mut new_tree, &new_root)?;

    Ok(new_tree)
}

pub fn get_subtree_with_leaves(
    tree: &Tree,
    leaves: &[usize],
    collapse_edges: bool,
) -> Result<Tree, TreeError> {
    let lca_id = get_lca(tree, leaves).ok_or(TreeError::GeneralError("Cannot get LCA"));
    let lca_id = lca_id?;
    let lca_node = tree.get(&lca_id)?.clone();

    let (dur, keep_set) = time(|| get_node_ids_for_leaf_set(tree, leaves, true));

    println!(
        "\t\tKeepset: {} for {} took {:?}",
        keep_set.len(),
        lca_id,
        dur
    );

    let mut new_tree = Tree::new();
    let mut new_root = Node::new();
    new_root.name = lca_node.name.clone();

    let (dur, res) =
        time(|| add_recursively_with_keep(tree, &mut new_tree, &keep_set, &lca_node, None));
    res?;

    println!("\t\tAdd recursively: {:?}", dur);

    if collapse_edges {
        let (dur, res) = time(|| tree_collapse_edges(&new_tree));
        new_tree = res?;
        println!("\t\tAdd recursively: {:?}", dur);
    }

    Ok(new_tree)
}

pub fn add_recursively(old_tree: &Tree, new_tree: &mut Tree, node: &Node) -> Result<(), TreeError> {
    new_tree.add(node.clone());

    for child_id in &node.children {
        let node = old_tree.get(&child_id)?;
        add_recursively(old_tree, new_tree, node)?;
    }

    Ok(())
}

pub fn add_recursively_with_keep(
    old_tree: &Tree,
    new_tree: &mut Tree,
    keep: &HashSet<usize>,
    old_node: &Node,
    new_parent_id: Option<&NodeId>,
) -> Result<(), TreeError> {
    let mut new_node = Node::new();
    new_node.name = old_node.name.clone();

    let new_node_id = match new_parent_id {
        Some(pid) => new_tree.add_child(new_node, *pid, old_node.parent_edge)?,
        None => new_tree.add(new_node),
    };

    // Stops if leaf-node
    for child_id in &old_node.children {
        // eprintln!("Is contained: {} ? {} (size: {})", child_id, keep.contains(child_id), keep.len());
        if !keep.contains(&child_id) {
            continue;
        };

        let child_node = old_tree.get(child_id)?;
        add_recursively_with_keep(old_tree, new_tree, keep, child_node, Some(&new_node_id))?;
    }

    Ok(())
}

pub fn get_lca(tree: &Tree, ids: &[usize]) -> Option<usize> {
    if ids.is_empty() {
        return None;
    }

    let mut lcas: HashSet<usize> = ids.iter().map(|x| *x).collect::<HashSet<usize>>();

    let mut tmp = HashSet::<usize>::new();
    if lcas.len() > 1 {
        while lcas.len() != 1 {
            lcas.iter()
                .tuple_combinations()
                .map(|(id1, id2)| tree.get_common_ancestor(id1, id2).unwrap())
                .collect_into(&mut tmp);
            swap(&mut lcas, &mut tmp);
            tmp.clear();
        }
    }
    Some(*lcas.iter().next().unwrap())
}

pub fn tree_collapse_edges_helper(
    new_tree: &mut Tree,
    active_node_id: &NodeId,
    tree: &Tree,
    current_id: &NodeId,
) -> Result<(), TreeError> {
    let children = (&tree.get(current_id)?).children.clone();

    if children.len() >= 1 {
        for child_id in children {
            let child_children = (&tree.get(&child_id)?).children.clone();

            if child_children.len() == 1 {
                // get first node that has more than one child
                // iterate while has only one child
                let mut end_id = child_children[0];

                let mut edge_length_sum = (&tree.get(&child_id)?).parent_edge.unwrap_or(0f64)
                    + (&tree.get(&end_id)?).parent_edge.unwrap_or(0f64);
                let mut intermediate_ids = vec![child_id];
                loop {
                    let current_children = (&tree.get(&end_id)?).children.clone();
                    if current_children.len() == 1 {
                        intermediate_ids.push(end_id);
                        end_id = current_children[0];
                        edge_length_sum += (&tree.get(&end_id)?).parent_edge.unwrap_or(0f64);
                        continue;
                    } else {
                        break;
                    }
                }

                let old_node = tree.get(&end_id)?;
                let mut new_node = Node::new();
                new_node.name = old_node.name.clone();
                let new_id =
                    new_tree.add_child(new_node, *active_node_id, Some(edge_length_sum))?;

                tree_collapse_edges_helper(new_tree, &new_id, tree, &end_id)?;
            } else {
                let old_node = tree.get(&child_id)?;
                let mut new_node = Node::new();
                new_node.name = old_node.name.clone();
                let new_id = new_tree.add_child(new_node, *active_node_id, old_node.parent_edge)?;
                tree_collapse_edges_helper(new_tree, &new_id, tree, &child_id)?;
            }
        }
    }

    Ok(())
}

pub fn tree_collapse_edges(tree: &Tree) -> Result<Tree, TreeError> {
    // To compress, the condition is that consecutive nodes have only one child each
    let mut new_tree = Tree::new();
    let root_node_id = tree.get_root()?;
    let root_node = tree.get(&root_node_id)?;

    let mut new_root = Node::new();
    new_root.name = root_node.name.clone();
    let new_root_id = new_tree.add(new_root);

    let _ = tree_collapse_edges_helper(&mut new_tree, &new_root_id, tree, &root_node_id);

    Ok(new_tree)
}

pub fn wrap_names(tree: &mut Tree) -> Result<(), TreeError> {
    let mut node_ids: Vec<_> = tree.get_nodes().clone();
    node_ids.extend(tree.get_leaves());

    for node_id in node_ids {
        // Get a mutable reference to the node
        if let Ok(node) = tree.get_mut(&node_id) {
            if let Some(name) = &mut node.name {
                *name = format!("\"{}\"", name); // Replace the original name
            }
        }
    }

    Ok(())
}

#[derive(Clone, Debug, PartialEq)]
pub struct NeighborDist {
    pub id: usize,
    pub name: String,
    pub distance: f64,
}
pub fn closest_neighbor<'a>(tree: &Tree, id: &NodeId) -> Result<NeighborDist, TreeError> {
    if tree.get_leaf_names().len() < 2 {
        return Err(TreeError::GeneralError("Tree has only 2 leaves"));
    }

    let dm = tree.distance_matrix()?;

    let node = tree.get_root()?;

    let name = tree
        .get(id)?
        .name
        .as_ref()
        .ok_or(TreeError::UnnamedLeaves)?;

    let closest = tree
        .get_leaf_names()
        .into_iter()
        .filter(|x| x.is_some() && x.as_ref().unwrap() != name)
        .map(|s| s.unwrap())
        .map(|leaf_name| (leaf_name.clone(), dm.get(&name, &leaf_name).unwrap()))
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();

    let node = tree.get_by_name(&closest.0).unwrap();

    Ok(NeighborDist {
        id: node.id,
        name: closest.0,
        distance: *closest.1,
    })
}

// pub fn closest_neighbor_str<'a>(tree: &Tree, id: &str) -> Result<(NodeId, f64), TreeError> {
//     let dm = tree.distance_matrix()?;

//     let node = tree.get_root()?;

//     let name = tree.get(id)?.name.as_ref().ok_or(TreeError::UnnamedLeaves)?;

//     println!("DM: {:?}", dm);

//     for leaf_name in tree.get_leaf_names() {
//         if let Some(leaf_name) = leaf_name {
//             if &leaf_name == name { continue };
//             println!("{} {} => {:?}", name, leaf_name, dm.get(&name, &leaf_name));
//         }
//     }

//     Ok((node, 0f64))
// }
