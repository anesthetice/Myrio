// Modules
mod misc;
mod run;
mod tree;

// Re-exports
use std::path::PathBuf;

// Imports
use anyhow::Context;
use anyhow::bail;
use clap::{Arg, ArgAction, ArgMatches, ColorChoice, Command, value_parser as vparser};
use console::style;
use indicatif::MultiProgress;
use indoc::printdoc;
use itertools::Itertools;
pub(in crate::app) use misc::{process_misc as misc, subcommand as misc_subcommand};
use myrio_core::{
    data::{MyrSeq, SFVec},
    tax::store::TaxTreeStore,
};
pub(in crate::app) use run::{process as run, subcommand as run_subcommand};
pub(in crate::app) use tree::{process_tree as tree, subcommand as tree_subcommand};

use crate::app::{App, config::Config};

fn gather_trees(
    mat: &mut ArgMatches,
    arg_id: &str,
) -> anyhow::Result<Vec<PathBuf>> {
    let mut tree_filepaths: Vec<PathBuf> = Vec::new();
    for pathbuf in mat.remove_many::<PathBuf>(arg_id).unwrap() {
        if pathbuf.is_dir() {
            pathbuf
                .read_dir()?
                .filter_map(Result::ok)
                .map(|dir_entry| dir_entry.path())
                .filter(|pathbuf| {
                    pathbuf.extension().and_then(|s| s.to_str()) == Some(TaxTreeStore::FILE_EXTENSION)
                })
                .for_each(|pathbuf| tree_filepaths.push(pathbuf));
        } else {
            if !pathbuf.is_file() {
                bail!("Invalid filepath: '{}'", pathbuf.display())
            }
            tree_filepaths.push(pathbuf);
        }
    }
    Ok(tree_filepaths)
}
