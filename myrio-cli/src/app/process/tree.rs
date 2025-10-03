use super::*;

pub fn subcommand() -> Command {
    let tree_new_subcommand = Command::new("new")
        .about("Create a new `.myrtree` database from a `.fasta` precursor file")
        .arg(
            Arg::new("gene")
                .help("The name of the gene")
                .required(true)
                .short('g')
                .long("gene")
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("input")
                .help("The `.fasta` file to use as input")
                .required(true)
                .short('i')
                .long("input")
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("output")
                .help("The output path/filepath to use")
                .required(false)
                .short('o')
                .long("output")
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("pre-compute")
                .help("The k-mer counts to pre-compute and store for one or more `k` ")
                .required(false)
                .num_args(1..101)
                .short('k')
                .long("pre-compute")
                .value_parser(vparser!(usize))
                .action(ArgAction::Set),
        );

    let tree_expand_subcommand = Command::new("expand")
        .about("Pre-compute k-mer counts for one or more `.myrtree` files")
        .arg(
            Arg::new("trees")
                .help("The one or more `.myrtree` reference databases")
                .long_help("The one or more `.myrtree` reference databases, also accepts directories")
                .required(true)
                .num_args(1..101)
                .short('t')
                .visible_short_alias('r')
                .long("trees")
                .visible_aliases(["refs", "references", "db"])
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("pre-compute")
                .help("The k-mer counts to pre-compute and store for one or more `k`")
                .long_help("The k-mer counts to pre-compute and store for one or more `k`, defaults to the config `search.k` if unset")
                .required(false)
                .num_args(1..101)
                .short('k')
                .long("pre-compute")
                .value_parser(vparser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("nb-resamples")
                .help("The number of times to resample/downsample each FASTA sequence")
                .long_help("The number of times to resample/downsample each FASTA sequence, defaults to `fasta_bootstrapping_nb_resamples_default` as specified in the config")
                .required(false)
                .short('n')
                .long("nb-resamples")
                .value_parser(vparser!(usize))
                .action(ArgAction::Set),

        );

    let tree_shrink_subcommand = Command::new("shrink")
        .about("Discard the pre-computed k-mer counts for one or more `.myrtree` files")
        .arg(
            Arg::new("trees")
                .help("The one or more `.myrtree` files to shrink")
                .long_help("The one or more `.myrtree` files to shrink (discard pre-computed k-mer counts), accepts filepaths as well as directories that will be searched for `.myrtree` files")
                .required(true)
                .num_args(1..101)
                .short('t')
                .visible_short_alias('i')
                .long("trees")
                .visible_aliases(["input"])
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set),
        );

    Command::new("tree").subcommands([tree_new_subcommand, tree_expand_subcommand, tree_shrink_subcommand])
}

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
