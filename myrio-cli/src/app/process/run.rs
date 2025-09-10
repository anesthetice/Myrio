use std::ops::{Div, DivAssign};

use indicatif::MultiProgress;
use itertools::Itertools;
use myrio_core::{
    clustering::ClusteringParameters,
    data::Float,
    similarity::{SimFunc, Similarity},
    tax::{
        compute::{CacheOptions, TaxTreeCompute},
        results::TaxTreeResults,
    },
};
use rand::SeedableRng;

use super::*;

pub fn process_run(
    mut mat: ArgMatches,
    config: &Config,
) -> anyhow::Result<()> {
    let input_filepath: PathBuf = mat.remove_one("input").unwrap();
    let tree_filepaths = gather_trees(&mut mat, "trees")?;

    let cluster_k = config.cluster.k;
    let cluster_t1 = config.cluster.t1;
    let cluster_t2 = config.cluster.t2;
    let cluster_similarity: Similarity = config.cluster.similarity.into();

    let fingerprint_nb_subsamples = config.nb_subsamples_for_tree_centroid;
    let fingerprint_fasta_nb_resamples = config.fasta_bootstrapping_nb_resamples_for_tree_centroid;
    let cache_opt = if mat.get_flag("cache-counts") {
        CacheOptions::Enabled { zstd_compression_level: config.zstd_compression_level }
    } else {
        CacheOptions::Disabled
    };

    let search_k = mat.remove_one::<usize>("k-search").unwrap_or(config.search.k);
    let search_t1 = config.search.t1;
    let search_similarity: Similarity = config.search.similarity.into();

    let myrseqs: Vec<MyrSeq> = myrio_core::io::read_fastq_from_file(input_filepath)?;

    // gather compute trees
    let mut rng = rand::rngs::StdRng::try_from_os_rng()?;
    let ttcompute_vec = tree_filepaths
        .iter()
        .map(|filepath| {
            let multi = MultiProgress::new();
            let spinner = myrio_core::utils::simple_spinner(
                Some(format!("Loading and processing '{}'", filepath.display())),
                Some(200),
                Some(&multi),
            );

            let ttstore = TaxTreeStore::load_from_file(filepath)?;

            let ttcompute = TaxTreeCompute::from_store_tree(
                ttstore,
                cluster_k,
                fingerprint_nb_subsamples,
                fingerprint_fasta_nb_resamples,
                search_k,
                config.fasta_bootstrapping_nb_resamples_default,
                cache_opt,
                &mut rng,
                Some(&multi),
            );
            spinner.finish();
            ttcompute
        })
        .collect::<Result<Vec<TaxTreeCompute>, myrio_core::tax::Error>>()?;

    let spinner = myrio_core::utils::simple_spinner(Some("Clustering...".to_string()), Some(200), None);

    #[rustfmt::skip]
    let nb_clusters = match mat.remove_one::<usize>("nb-clusters") {
        Some(nb_clusters) => {
            if nb_clusters < tree_filepaths.len() { bail!("Expected `nb-clusters` to be greater or equal to the number of specified `.myrtree` files") }
            nb_clusters
        }
        None => tree_filepaths.len(),
    };

    let intial_centroids = ttcompute_vec
        .iter()
        .map(TaxTreeCompute::get_kmer_normalized_counts_fingerprint)
        .cloned()
        .collect_vec();

    let cluster_params = ClusteringParameters {
        k: cluster_k,
        t1: cluster_t1,
        t2: cluster_t2,
        similarity: cluster_similarity,
        intial_centroids,
        expected_nb_of_clusters: nb_clusters,
        eta_improvement: config.cluster.eta_improvement,
        nb_iters_max: config.cluster.nb_iters_max,
        silhouette_trimming: config.cluster.silhouette_trimming,
    };

    let queries = myrio_core::clustering::cluster(myrseqs, cluster_params)
        .clusters
        .into_iter()
        .map(|myrseqs| {
            let elements = myrseqs
                .into_iter()
                .map(|myrseq| myrseq.compute_kmer_counts(search_k, search_t1).0)
                .collect_vec();
            myrio_core::clustering::compute_cluster_centroid(&elements)
        })
        .collect_vec();

    spinner.finish_with_message("Finished clustering");

    for (ttcompute, query) in ttcompute_vec.into_iter().zip(queries) {
        let ttresults_full = TaxTreeResults::from_compute_tree(
            query,
            ttcompute,
            search_similarity,
            config.search.lambda_leaf,
            config.search.lambda_branch,
            config.search.mu,
            config.search.gamma,
            config.search.delta,
            config.search.epsilon,
            None,
        )?;

        /*
        use std::io::Write;
        let mut file =
            std::fs::OpenOptions::new().write(true).create(true).truncate(true).open("./output.txt")?;
        writeln!(file, "{ttresults_full}")?;
        file.sync_all()?;
        drop(file);
        break;
        */

        let ttresults_best = ttresults_full.cut(config.search.nb_best_display);

        println!("{}", ttresults_best.core);
    }

    Ok(())
}
