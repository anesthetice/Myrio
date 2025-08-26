use std::path::PathBuf;

use clap::ArgMatches;
use myrio_core::tax::store::TaxTreeStore;

use crate::app::config::Config;

pub fn process_tree(
    mut mat: ArgMatches,
    config: &Config,
) -> anyhow::Result<()> {
    let Some((subcommand, sub_mat)) = mat.remove_subcommand() else {
        return Ok(());
    };

    match subcommand.as_str() {
        "new" => process_tree_new(sub_mat, config),
        _ => Ok(()),
    }
}

fn process_tree_new(
    mut mat: ArgMatches,
    config: &Config,
) -> anyhow::Result<()> {
    let gene: String = mat.remove_one("gene").unwrap();
    let input: PathBuf = mat.remove_one("input").unwrap();

    let output = mat.remove_one::<PathBuf>("output").map_or(Ok(None), |mut pathbuf| {
        if pathbuf.is_dir() {
            let input_filename = input.file_name().ok_or_else(|| {
                anyhow::anyhow!(
                    "Provided input filepath does not contain a valid filename, '{}'",
                    input.display()
                )
            })?;
            pathbuf.set_file_name(input_filename);
            pathbuf.set_extension(TaxTreeStore::FILE_EXTENSION);
        }
        anyhow::Ok(Some(pathbuf))
    })?;

    Ok(())
}
