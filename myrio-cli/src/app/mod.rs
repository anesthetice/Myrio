// Modules
pub mod cli;
pub mod config;
mod process;

// Imports
use anyhow::Context;
use clap::ColorChoice;
use config::Config;
use directories::ProjectDirs;

pub struct App {
    config: Config,
}

impl App {
    #[cfg(not(debug_assertions))]
    const NAME: &'static str = "myrio";
    #[cfg(debug_assertions)]
    const NAME: &'static str = "myrio-dev";

    const VERSION: &'static str = env!("CARGO_PKG_VERSION");

    pub fn load() -> anyhow::Result<Self> {
        let dirs = ProjectDirs::from("", "", Self::NAME)
            .ok_or_else(|| anyhow::anyhow!("Failed to get project directories"))?;

        let _conf_dir = dirs.config_dir();
        if !_conf_dir.exists() {
            std::fs::create_dir_all(_conf_dir).with_context(|| {
                format!("Failed to create the config directory at `{}`", _conf_dir.display())
            })?;
        }

        let config = config::Config::load(&_conf_dir.join(Config::FILENAME));

        Ok(Self { config })
    }

    pub fn run(self) -> anyhow::Result<()> {
        let mut mat = cli::build_cli().get_matches();

        if mat.get_flag("version") {
            println!("myrio {}", Self::VERSION);
            return Ok(());
        }

        let color_choice = mat.remove_one::<String>("color").unwrap();
        let color_choice = match color_choice.as_str() {
            "always" => {
                console::set_colors_enabled(true);
                console::set_colors_enabled_stderr(true);
                ColorChoice::Always
            }
            "auto" => ColorChoice::Auto,
            "never" => {
                console::set_colors_enabled(false);
                console::set_colors_enabled_stderr(false);
                ColorChoice::Never
            }
            _ => unreachable!(),
        };

        let Some((subcommand, sub_mat)) = mat.remove_subcommand() else {
            return Ok(());
        };

        match subcommand.as_str() {
            "run" => process::run(sub_mat, &self.config, color_choice),
            "tree" => process::tree(sub_mat, &self.config),
            "misc" => process::misc(sub_mat, &self.config),
            _ => Ok(()),
        }
    }
}
