use std::path::PathBuf;

use clap::ArgMatches;
use myrio_core::tax::store::TaxTreeStore;

use crate::app::App;

pub fn process_tree(
    mut mat: ArgMatches,
    app_ref: &App,
) -> anyhow::Result<()> {
    let Some((subcommand, sub_mat)) = mat.remove_subcommand() else {
        return Ok(());
    };

    match subcommand.as_str() {
        "new" => process_tree_new(sub_mat, app_ref),
        _ => Ok(()),
    }
}

fn process_tree_new(
    mut mat: ArgMatches,
    app_ref: &App,
) -> anyhow::Result<()> {
    let gene: String = mat.remove_one("gene").unwrap();
    let input: PathBuf = mat.remove_one("input").unwrap();

    #[rustfmt::skip]
    let output = mat
        .remove_one::<PathBuf>("output")
        .map_or(Ok(None), |mut pathbuf| {
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

    let pre_compute_k_values: Option<Vec<usize>> = mat.remove_many("pre-compute").map(|vals| vals.collect());
    let pre_compute_kcounts = pre_compute_k_values
        .as_deref()
        .map(|k_values| (k_values, app_ref.config.fasta_expansion_max_consecutive_N_before_gap));

    let tree = TaxTreeStore::load_from_fasta_file(input, output, gene, pre_compute_kcounts)?;

    tree.encode_to_file(
        app_ref.config.zstd_compression_level,
        app_ref.config.zstd_multithreading_flag,
        app_ref.nb_threads_available,
    )?;

    Ok(())
}
