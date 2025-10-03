use super::*;

pub fn subcommand() -> Command {
    let misc_generate_shell_completions_subcommand = Command::new("generate-shell-completions")
        .about("Generate completions for your desired shell")
        .long_about(
            "This subcommand is used to generate shell completions for the selected shell, outputs to stdout",
        )
        .arg(
            Arg::new("shell")
                .index(1)
                .required(true)
                .help("The shell to target")
                .action(ArgAction::Set)
                .value_parser(vparser!(clap_complete::Shell)),
        );

    let misc_get_config_location =
        Command::new("get-config-location").about("Display the filepath of the configuration");

    Command::new("misc").subcommands([misc_generate_shell_completions_subcommand, misc_get_config_location])
}

pub fn process_misc(
    mut mat: ArgMatches,
    _config: &Config,
) -> anyhow::Result<()> {
    let Some((subcommand, sub_mat)) = mat.remove_subcommand() else {
        return Ok(());
    };

    match subcommand.as_str() {
        "generate-shell-completions" => process_misc_generate_shell_completions(sub_mat, _config),
        "get-config-location" => {
            let conf_filepath = directories::ProjectDirs::from("", "", App::NAME)
                .ok_or_else(|| anyhow::anyhow!("Failed to get project directories"))?
                .config_dir()
                .join(Config::FILENAME);

            println!("Configuration filepath: {}", conf_filepath.display());
            Ok(())
        }
        _ => Ok(()),
    }
}

fn process_misc_generate_shell_completions(
    mut mat: ArgMatches,
    _config: &Config,
) -> anyhow::Result<()> {
    if let Some(generator) = mat.remove_one::<clap_complete::Shell>("shell") {
        let mut cmd = crate::app::cli::build_cli();
        let bin_name = cmd.get_name().to_string();
        clap_complete::generate(generator, &mut cmd, bin_name, &mut std::io::stdout());
        Ok(())
    } else {
        bail!("No valid shell (generator) provided")
    }
}
