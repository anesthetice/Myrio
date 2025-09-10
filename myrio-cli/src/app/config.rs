use std::io::{Read, Write};

use myrio_core::data::Float;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub zstd_compression_level: i32,
    pub nb_subsamples_for_tree_centroid: usize,
    pub fasta_bootstrapping_nb_resamples_default: usize,
    pub fasta_bootstrapping_nb_resamples_for_tree_centroid: usize,
    pub cluster: ClusterConfig,
    pub search: SearchConfig,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            zstd_compression_level: 9,
            nb_subsamples_for_tree_centroid: 512,
            fasta_bootstrapping_nb_resamples_default: 32,
            fasta_bootstrapping_nb_resamples_for_tree_centroid: 32,
            cluster: ClusterConfig::default(),
            search: SearchConfig::default(),
        }
    }
}

impl Config {
    pub const DEFAULT_WITH_COMMENTS: &'static str = indoc::indoc! {r#"
        # The compression level used by zstd, levels currently range from 1 to 22
        zstd_compression_level = 6

        # The number of subsamples (i.e., leaves) used to compute the kmer counts centroid of a tree
        nb_subsamples_for_tree_centroid = 512

        # The number of times resampling (techincally downsampling) is performed on a fasta sequence when computing kmer counts by default
        fasta_bootstrapping_nb_resamples_default = 32

        # The number of times resampling (techincally downsampling) is performed on a fasta sequence when computing kmer counts for the tree centroid
        fasta_bootstrapping_nb_resamples_for_tree_centroid = 32

        [cluster]
        k = 6 # k-mer size to use in clustering (i.e., value of `k` itself)
        t1 = 0.2 # excludes k-mers from the count for which the product of their per-base probability of being correctly identified is lower than this value`
        t2 = 20 # reject sequences that have fewer than `t2` valid k-mers
        similarity = "Cosine" # similarity function to use, available: ["Cosine", "JacardTanimoto", "Overlap"]
        eta_improvement = 0.0001 # fraction used to decide when to stop the main clustering step
        nb_iters_max = 20 # maximum number of iterations the main clustering step can usec
        # E.g., `silhouette_trimming = 1.4` means silhouette trimming is enabled, sequences with a lower silhouette score than `mean - std * 1.4` are discarded
        # To disabled silhouette trimming, comment out the line with `#` (improves clustering performance by a lot)
        silhouette_trimming = 1.4

        [search]
        k = 18 # k-mer size to use in searching (i.e., value of `k` itself)
        t1 = 0.1 # see `cluster.t1` above
        similarity = "Cosine" # see `cluster.similarity` above
        nb_best_analysis = 100
        lambda_leaf = 2.0
        lambda_branch = 20.0
        mu = 1.1
        gamma = 2.5
        delta = 0.9
        epsilon = 0.55
        nb_best_display = 7
    "#};

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

/// From `myrio-core/src/similarity.rs`
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Similarity {
    Cosine,
    JacardTanimoto,
    Overlap,
}

#[allow(clippy::from_over_into)]
impl Into<myrio_core::similarity::Similarity> for Similarity {
    fn into(self) -> myrio_core::similarity::Similarity {
        match self {
            Self::Cosine => myrio_core::similarity::Similarity::Cosine,
            Self::JacardTanimoto => myrio_core::similarity::Similarity::JacardTanimoto,
            Self::Overlap => myrio_core::similarity::Similarity::Overlap,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename = "clustering")]
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
#[serde(rename = "searching")]
pub struct SearchConfig {
    pub k: usize,
    pub t1: Float,
    pub similarity: Similarity,
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
