// Imports
use clap::{Arg, ArgAction, Command, builder::Styles};

use crate::app::process;

#[rustfmt::skip]
pub fn build_cli() -> Command {
    Command::new("myrio")
        .color(clap::ColorChoice::Auto)
        .styles(Styles::styled())
        .arg(
            Arg::new("version")
                .required(false)
                .short('v')
                .long("version")
                .action(ArgAction::SetTrue)
        )
        .arg(
            Arg::new("color")
                .required(false)
                .long("color")
                .value_parser(["always", "auto", "never"])
                .default_value("auto")
                .action(ArgAction::Set)
        )
        .subcommands([
            process::run_subcommand(),
            process::tree_subcommand(),
            process::misc_subcommand(),
        ])
}
