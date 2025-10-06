// Imports
use std::io::{Read, Write};

use myrio_core::{data::Float, similarity::Similarity, tax::clade::Rank};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub zstd_compression_level: i32,
    pub fastq_min_length: usize,
    pub fastq_min_mean_qual: Float,
    pub fastq_max_qual: u8,
    pub fasta_nb_resamples: usize,

    pub cluster: ClusterConfig,
    pub fingerprint: FingerprintConfig,
    pub search: SearchConfig,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            zstd_compression_level: 6,
            fastq_min_length: 50,
            fastq_min_mean_qual: 12.0,
            fastq_max_qual: 255,
            fasta_nb_resamples: 32,
            cluster: ClusterConfig::default(),
            fingerprint: FingerprintConfig::default(),
            search: SearchConfig::default(),
        }
    }
}

impl Config {
    pub const DEFAULT_WITH_COMMENTS: &'static str = include_str!("config_default.toml");
    pub const FILENAME: &str = "myrio.conf.toml";

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
                    //match config.to_file(filepath) {
                    match Self::default_with_comments_to_file(filepath) {
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

    #[allow(unused)]
    fn to_file(
        &self,
        filepath: &std::path::Path,
    ) -> anyhow::Result<()> {
        let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(filepath)?;
        file.write_all(toml::to_string_pretty(&self)?.as_bytes())?;
        file.sync_all()?;
        Ok(())
    }

    fn default_with_comments_to_file(filepath: &std::path::Path) -> anyhow::Result<()> {
        let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(filepath)?;
        file.write_all(Self::DEFAULT_WITH_COMMENTS.as_bytes())?;
        file.sync_all()?;
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClusterConfig {
    pub k: usize,
    pub t1: Float,
    pub t2: usize,
    pub similarity: Similarity,
    pub eta_improvement: Float,
    pub nb_iters_max: usize,
    pub silhouette_trimming: Option<Float>,
}

impl Default for ClusterConfig {
    fn default() -> Self {
        Self {
            k: 6,
            t1: 0.2,
            t2: 20,
            similarity: Similarity::Cosine,
            eta_improvement: 1E-4,
            nb_iters_max: 20,
            silhouette_trimming: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FingerprintConfig {
    pub local_rank: Rank,
    pub similarity: Similarity,
    pub local_nb_subsamples: usize,
    pub local_fasta_nb_resamples: usize,
    pub alpha: Float,
}

impl Default for FingerprintConfig {
    fn default() -> Self {
        Self {
            local_rank: Rank::Phylum,
            similarity: Similarity::Cosine,
            local_nb_subsamples: 128,
            local_fasta_nb_resamples: 32,
            alpha: 2.5,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchConfig {
    pub k: usize,
    pub t1: Float,
    pub similarity: Similarity,
    pub max_amount_of_leaves_per_branch: usize,
    pub nb_best_analysis: usize,
    pub lambda_leaf: Float,
    pub lambda_branch: Float,
    pub mu: Float,
    pub gamma: Float,
    pub delta: Float,
    pub epsilon: Float,
    pub nb_best_display: usize,
}

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            k: 18,
            t1: 0.1,
            similarity: Similarity::Cosine,
            max_amount_of_leaves_per_branch: 15,
            nb_best_analysis: 100,
            lambda_leaf: 2.0,
            lambda_branch: 20.0,
            mu: 1.1,
            gamma: 2.5,
            delta: 0.9,
            epsilon: 0.55,
            nb_best_display: 7,
        }
    }
}
