use std::io::{Read, Write};

use serde::{Deserialize, Serialize};

#[allow(non_snake_case)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub zstd_compression_level: i32,
    pub nb_bootstrap_resamples: usize,
    pub cluster_k_default: usize,
    pub search_k_default: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            zstd_compression_level: 9,
            nb_bootstrap_resamples: 32,
            cluster_k_default: 6,
            search_k_default: 18,
        }
    }
}

impl Config {
    pub fn load(filepath: &std::path::Path) -> Self {
        match Self::from_file(filepath) {
            Ok(config) => config,
            Err(err) => {
                eprintln!("Warning: failed to load configuration from file, '{err}'");
                let config = Config::default();
                let Ok(downcast_error) = err.downcast::<std::io::Error>() else {
                    return config;
                };
                if downcast_error.kind() == std::io::ErrorKind::NotFound {
                    match config.to_file(filepath) {
                        Ok(()) => eprintln!(
                            "Warning: Created default configuration file, at '{}'",
                            filepath.display()
                        ),
                        Err(error) => eprintln!(
                            "Warning: Failed to create default configuration file, at '{}', caused by '{}'",
                            filepath.display(),
                            error
                        ),
                    }
                }
                config
            }
        }
    }

    fn from_file(filepath: &std::path::Path) -> anyhow::Result<Self> {
        let mut buffer: Vec<u8> = Vec::new();
        std::fs::OpenOptions::new()
            .create(false)
            .read(true)
            .open(filepath)?
            .read_to_end(&mut buffer)?;
        Ok(toml::from_slice(&buffer)?)
    }

    fn to_file(
        &self,
        filepath: &std::path::Path,
    ) -> anyhow::Result<()> {
        let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(filepath)?;
        file.write_all(toml::to_string_pretty(&self)?.as_bytes())?;
        file.sync_all()?;
        Ok(())
    }
}
