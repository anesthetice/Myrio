use super::*;

pub fn process_tree(
    mut mat: ArgMatches,
    config: &Config,
) -> anyhow::Result<()> {
    let Some((subcommand, sub_mat)) = mat.remove_subcommand() else {
        return Ok(());
    };

    match subcommand.as_str() {
        "new" => process_tree_new(sub_mat, config),
        "expand" => process_tree_expand(sub_mat, config),
        "shrink" => process_tree_shrink(sub_mat, config),
        _ => Ok(()),
    }
}

fn process_tree_new(
    mut mat: ArgMatches,
    config: &Config,
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

    let mut ttstore = TaxTreeStore::load_from_fasta_file(input, output, gene)?;

    if let Some(k_values_to_precompute) = mat.remove_many::<usize>("pre-compute") {
        for k in k_values_to_precompute.unique() {
            ttstore.compute_and_append_kmer_counts(k, config.fasta_nb_resamples, None);
        }
    }

    ttstore.save_to_file(config.zstd_compression_level, None)?;

    Ok(())
}

fn process_tree_expand(
    mut mat: ArgMatches,
    config: &Config,
) -> anyhow::Result<()> {
    let tree_filepaths = gather_trees(&mut mat, "trees")?;

    let k_values_to_precompute_vec = mat
        .remove_many::<usize>("pre-compute")
        .map(|vals| vals.collect_vec()) // no need for unique here
        .unwrap_or_else(|| vec![config.search.k]);

    let nb_resamples = mat.remove_one::<usize>("nb-resamples").unwrap_or(config.fasta_nb_resamples);

    for filepath in tree_filepaths.into_iter() {
        let spinner = myrio_core::utils::simple_spinner(
            Some(format!("Loading '{}'", filepath.display())),
            Some(200),
            None,
        );
        let mut ttstore = TaxTreeStore::load_from_file(filepath)?;
        spinner.finish();

        for &k in k_values_to_precompute_vec.iter() {
            if !ttstore.k_precomputed.contains(&k) {
                ttstore.compute_and_append_kmer_counts(k, nb_resamples, None);
            }
        }

        ttstore.save_to_file(config.zstd_compression_level, None)?;
    }
    Ok(())
}

fn process_tree_shrink(
    mut mat: ArgMatches,
    config: &Config,
) -> anyhow::Result<()> {
    let tree_filepaths = gather_trees(&mut mat, "trees")?;

    for filepath in tree_filepaths.into_iter() {
        let spinner = myrio_core::utils::simple_spinner(
            Some(format!("Loading '{}'", filepath.display())),
            Some(200),
            None,
        );
        let mut ttstore = TaxTreeStore::load_from_file(filepath)?;
        spinner.finish();

        ttstore.shrink();
        ttstore.save_to_file(config.zstd_compression_level, None)?;
    }
    Ok(())
}
