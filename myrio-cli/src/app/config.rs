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
    pub const FILENAME: &str = "myrio.conf.toml";

    pub const DEFAULT_WITH_COMMENTS: &'static str = indoc::indoc! {r#"
        # The compression level used by zstd, levels currently range from 1 to 22
        zstd_compression_level = 6

        # The minimum length a fastq record needs to not be discarded
        fastq_min_length = 50
        # The minimum mean quality (Phred score, typically between 0-40) a fastq record needs to not be discarded
        fastq_min_mean_qual = 12.0
        # The maximum quality score (Phred score) a fastq record can have before being discarded due to being anomalous
        fastq_max_qual = 255 # setting this value to 255 is basically the same as disabling it

        # The number of times resampling is performed on a fasta sequence when computing k-mer counts
        fasta_nb_resamples = 32


        [cluster]
        # The k-mer size to use in clustering (i.e., value of `k` itself)
        k = 6
        # Value that serves to exclude k-mers from being counted when `∏ P(correctly called) < t1`
        t1 = 0.2
        # Value that serves to exclude sequences that end up having fewer than `t2` valid k-mers
        t2 = 20
        # The similarity function to use, available: ["Cosine", "JacardTanimoto", "Overlap"]
        similarity = "Cosine"
        # A fraction used to decide when to stop the main clustering step
        eta_improvement = 0.0001
        # Maximum number of iterations the main clustering step can perform
        nb_iters_max = 20
        # E.g., `silhouette_trimming = 1.4` means silhouette trimming is enabled, sequences with a lower silhouette score than `mean - std * 1.4` are discarded
        # To disabled silhouette trimming, comment out the line with `#` (improves clustering performance by a lot)
        silhouette_trimming = 1.4


        [fingerprint]
        # Rank for which the local fingerprints will be computed, available: ["Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
        local_rank = "Class"
        # Same explanation as `cluster.similarity` seen above
        similarity = "Cosine"
        # The number of subsamples (i.e., leaves) used to compute the local fingerprint (representative k-mer counts) of a tree
        local_nb_subsamples = 128
        # The number of times resampling is performed on a fasta sequence when computing k-mer counts
        local_fasta_nb_resamples = 32
        # Mathematical parameter used to prioritize local fingerprints with a higher similarity score
        # equation: `local_fingerprint * exp(1 + α ⋅ score)`
        alpha = 2.5


        [search]
        # The k-mer size to use in searching (i.e., value of `k` itself)
        k = 18
        # Same explanation as `cluster.t1` seen above
        t1 = 0.1
        # Same explanation as `cluster.similarity` seen above
        similarity = "Cosine"
        # The maximum number of leaves per branch, if a branch contains more leaves than this value, then only the best are kept, the rest being discarded (from the pooling score computation but not the previous more basic ones)
        max_amount_of_leaves_per_branch = 15
        # The number of leaves to keep for analysis
        nb_best_analysis = 100

        # Mathematical parameters for softmax pooling (and small `n` penalty for leaves)
        # see `https://github.com/anesthetice/Myrio/blob/main/myrio-core/README.md` for more information
        lambda_leaf = 2.0
        lambda_branch = 20.0
        mu = 1.1

        # Mathematical parameters for confidence computation
        # see `https://github.com/anesthetice/Myrio/blob/main/myrio-core/README.md` for more information
        gamma = 2.5
        delta = 0.9
        epsilon = 0.55

        # Decides how many of the best leaves are displayed at the end
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
