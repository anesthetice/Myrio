use std::io::{Read, Write};

use serde::{Deserialize, Serialize};

#[allow(non_snake_case)]
#[derive(Debug, Clone, Deserialize)]
#[serde(from = "ConfigPrecursor")]
pub struct Config {
    pub zstd_compression_level: i32,
    pub zstd_multithreading_opt: Option<u32>,
    pub fasta_expansion_max_consecutive_N_before_gap: usize,
}

impl Config {
    pub fn load(filepath: &std::path::Path) -> Self {
        match Self::from_file(filepath) {
            Ok(config) => config,
            Err(err) => {
                eprintln!("Warning: failed to load configuration from file, '{err}'");
                let config = ConfigPrecursor::default();
                let Ok(downcast_error) = err.downcast::<std::io::Error>() else {
                    return config.into();
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
                config.into()
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
        Ok(ijson::from_value(&serde_json::from_slice(&buffer)?)?)
    }
}

impl From<ConfigPrecursor> for Config {
    fn from(value: ConfigPrecursor) -> Self {
        let zstd_multithreading_opt = if value.zstd_multithreading_flag {
            match std::thread::available_parallelism() {
                Ok(num) => Some(usize::from(num) as u32),
                Err(e) => {
                    eprintln!("Failed to get available parallelism for zstd, {e}");
                    None
                }
            }
        } else {
            None
        };
        eprintln!("{zstd_multithreading_opt:?}");
        Self {
            zstd_compression_level: value.zstd_compression_level,
            zstd_multithreading_opt,
            fasta_expansion_max_consecutive_N_before_gap: value.fasta_expansion_max_consecutive_N_before_gap,
        }
    }
}

#[allow(non_snake_case)]
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, rename = "config")]
pub struct ConfigPrecursor {
    pub zstd_compression_level: i32,
    pub zstd_multithreading_flag: bool,
    pub fasta_expansion_max_consecutive_N_before_gap: usize,
}

impl Default for ConfigPrecursor {
    fn default() -> Self {
        Self {
            zstd_compression_level: 15,
            zstd_multithreading_flag: true,
            fasta_expansion_max_consecutive_N_before_gap: 3,
        }
    }
}

impl ConfigPrecursor {
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
