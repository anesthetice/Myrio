use std::io::{Read, Write};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct Config {
    zstd_compression_level: i32,
}

impl Default for Config {
    fn default() -> Self {
        Self { zstd_compression_level: 17 }
    }
}

impl Config {
    pub fn load(filepath: &std::path::Path) -> Self {
        match Self::from_file(filepath) {
            Ok(config) => config,
            Err(err) => {
                eprintln!("Warning: Failed to load configuration from file, '{err}'");
                let config = Self::default();
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
        std::fs::OpenOptions::new().create(false).read(true).open(filepath)?.read_to_end(&mut buffer)?;
        Ok(ijson::from_value(&serde_json::from_slice(&buffer)?)?)
    }

    fn to_file(
        &self,
        filepath: &std::path::Path,
    ) -> anyhow::Result<()> {
        let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(filepath)?;

        file.write_all(&serde_json::to_vec_pretty(&ijson::to_value(self)?)?)?;
        file.flush()?;
        Ok(())
    }
}
