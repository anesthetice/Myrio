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
        ).arg(
            Arg::new("nb_clusters")
                .help("The number of clusters to expect")
                .long_help("The number of clusters to expect, defaults to the number of `.myrtree` files found")
                .required(false)
                .short('n')
                .long("nb-clusters")
                .value_parser(vparser!(usize))
                .action(ArgAction::Set)
        );

    let tree_new_subcommand = Command::new("new")
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

    let tree_shrink_subcommand = Command::new("shrink")
        .arg(
            Arg::new("trees")
                .help("The one or more `.myrtree` files to shrink")
                .long_help("The one or more `.myrtree` files to shrink (discard pre-computed kmer-freqs), accepts filepaths as well as directories that will be searched for `.myrtree` files")
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
            tree_shrink_subcommand,
        ]);



    Command::new("myrio").subcommands([
        run_subcommand,
        tree_subcommand,
    ])
}
