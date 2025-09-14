// Imports
use std::path::PathBuf;

use clap::{Arg, ArgAction, Command, value_parser as vparser};

#[rustfmt::skip]
pub fn build_cli() -> Command {
    let run_subcommand = Command::new("run")
        .arg(
            Arg::new("input")
                .help("The `.fastq` file to use as input")
                .required(true)
                //.num_args(1..101)
                .short('i')
                .visible_short_alias('q')
                .long("input")
                .visible_alias("query")
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set),
        )
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
            Arg::new("k-search")
                .help("The length of each k-mer (i.e., `k` itself) used for sequence comparison")
                .required(false)
                .short('k')
                .long("k-search")
                .value_parser(vparser!(usize))
                .action(ArgAction::Set)
        )
        .arg(
            Arg::new("output-csv")
                .help("Write results to a `.csv` file")
                .required(false)
                .long("output-csv")
                .visible_aliases(["csv", "csv-output"])
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set)
        )
        .arg(
            Arg::new("output-txt")
                .help("Write results to a `.txt` file")
                .required(false)
                .long("output-txt")
                .visible_aliases(["txt", "txt-output"])
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set)
        )
        .arg(
            Arg::new("cache-counts")
                .help("Flag that decides if newly-computed kmer counts are then cached")
                .required(false)
                .long("cache-counts")
                .action(ArgAction::SetTrue)
        )
        .arg(
            Arg::new("nb-clusters")
                .help("The number of clusters to expect")
                .long_help("The number of clusters to expect, defaults to the number of `.myrtree` files found")
                .required(false)
                .short('n')
                .long("nb-clusters")
                .value_parser(vparser!(usize))
                .action(ArgAction::Set)
        );
        /*
        .arg(
            Arg::new("no-initial-centroids")
                .help("Flag that prevents the use of intial centroids derived from gene trees when clustering, instead matching afterwards")
                .required(false)
                .long("no-initial-centroids")
                .action(ArgAction::SetTrue)
        );
        */

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

    let tree_subcommand = Command::new("tree")
        .subcommands([
            tree_new_subcommand,
            tree_expand_subcommand,
            tree_shrink_subcommand,
        ]);


    Command::new("myrio").subcommands([
        run_subcommand,
        tree_subcommand,
    ])
}
