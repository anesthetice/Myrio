use clap::ArgMatches;

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

    Ok(())
}
