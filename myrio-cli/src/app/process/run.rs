use std::{
    fs::{File, OpenOptions},
    io::Write,
    time::SystemTime,
};

use anyhow::Context;
use clap::ColorChoice;
use indicatif::MultiProgress;
use itertools::Itertools;
use myrio_core::{
    clustering::ClusteringParameters,
    similarity::{SimFunc, Similarity},
    tax::{
        compute::{CacheOptions, TaxTreeCompute},
        results::TaxTreeResults,
    },
};

use super::*;

pub fn process_run(
    mut mat: ArgMatches,
    config: &Config,
    color_choice: ColorChoice,
) -> anyhow::Result<()> {
    let input_filepath: PathBuf = mat.remove_one("input").unwrap();
    let tree_filepaths = gather_trees(&mut mat, "trees")?;

    let cluster_k = config.cluster.k;
    let cluster_t1 = config.cluster.t1;
    let cluster_t2 = config.cluster.t2;
    let cluster_similarity = config.cluster.similarity;

    let fingerprint_local_rank = config.fingerprint.local_rank;
    let fingerprint_similarity = config.fingerprint.similarity;
    let fingerprint_local_nb_subsamples = config.fingerprint.local_nb_subsamples;
    let fingerprint_local_fasta_nb_resamples = config.fingerprint.local_fasta_nb_resamples;
    let fingerprint_alpha = config.fingerprint.alpha;

    let cache_opt = if mat.get_flag("cache-counts") {
        CacheOptions::Enabled { zstd_compression_level: config.zstd_compression_level }
    } else {
        CacheOptions::Disabled
    };

    let search_k = mat.remove_one::<usize>("k-search").unwrap_or(config.search.k);
    let search_t1 = config.search.t1;
    let search_similarity: Similarity = config.search.similarity;

    // gather myrseqs and filter out the poor quality ones
    let myrseqs: Vec<MyrSeq> = {
        let mut output = myrio_core::io::read_fastq_from_file(input_filepath)?;
        let initial_count = output.len();
        output = MyrSeq::pre_process(
            output,
            config.fastq_min_length,
            config.fastq_min_mean_qual,
            config.fastq_max_qual,
        );
        let final_count = output.len();
        println!("Pre-processed {initial_count} records → discarded {}", initial_count - final_count);
        output
    };

    // gather compute trees
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
                fingerprint_local_rank,
                fingerprint_local_nb_subsamples,
                fingerprint_local_fasta_nb_resamples,
                config.fasta_nb_resamples,
                search_k,
                cache_opt,
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

    let queries: Vec<SFVec> = {
        spinner.set_message("Clustering from global fingerprints...");

        let initial_centroids_from_global_fingerprints =
            ttcompute_vec.iter().map(|ttcompute| &ttcompute.global_fingerprint).cloned().collect_vec();

        let cluster_params_for_global = ClusteringParameters {
            k: cluster_k,
            t1: cluster_t1,
            t2: cluster_t2,
            similarity: cluster_similarity,
            intial_centroids: initial_centroids_from_global_fingerprints,
            expected_nb_of_clusters: nb_clusters,
            eta_improvement: config.cluster.eta_improvement,
            nb_iters_max: config.cluster.nb_iters_max,
            silhouette_trimming: None, // No silhouette trimming yet no matter what is set in the config
        };

        let naive_queries = myrio_core::clustering::cluster(myrseqs.clone(), cluster_params_for_global)
            .clusters
            .into_iter()
            .map(|myrseqs| {
                let elements = myrseqs
                    .into_iter()
                    .map(|myrseq| myrseq.compute_kmer_counts(cluster_k, cluster_t1).0)
                    .collect_vec();
                myrio_core::clustering::compute_cluster_centroid(&elements)
            })
            .collect_vec();

        spinner.set_message("Re-clustering from local fingerprints...");
        let simfunc: SimFunc = fingerprint_similarity.to_simfunc(true);
        let initial_centroids_from_local_fingerprints = naive_queries
            .into_iter()
            .zip(ttcompute_vec.iter())
            .map(|(naive_query, ttcompute)| {
                //println!();
                ttcompute
                    .local_fingerprints
                    .iter()
                    .map(|local_fingerprint| {
                        let score = *simfunc(local_fingerprint, &naive_query);
                        let fingerprint_adjustor = (1.0 + fingerprint_alpha * score).exp();
                        // println!("{score:.2} → {fingerprint_adjustor:.2}");
                        local_fingerprint * fingerprint_adjustor
                    })
                    .sum::<SFVec>()
                    .into_normalized_l2()
            })
            .collect_vec();

        let cluster_params_for_local = ClusteringParameters {
            k: cluster_k,
            t1: cluster_t1,
            t2: cluster_t2,
            similarity: cluster_similarity,
            intial_centroids: initial_centroids_from_local_fingerprints,
            expected_nb_of_clusters: nb_clusters,
            eta_improvement: config.cluster.eta_improvement,
            nb_iters_max: config.cluster.nb_iters_max,
            silhouette_trimming: config.cluster.silhouette_trimming,
        };

        let clusters = myrio_core::clustering::cluster(myrseqs, cluster_params_for_local).clusters;

        if let Some(pathbuf) = mat.remove_one::<PathBuf>("save-clusters") {
            if !pathbuf.is_dir() {
                std::fs::create_dir(&pathbuf).with_context(|| {
                    format!("Failed to create `{}`, please make sure the parent path exists or the right permissions are granted.", pathbuf.display())
                })?;
            }

            let mut ttcompute_iter = ttcompute_vec.iter();
            for (idx, cluster) in clusters.iter().enumerate() {
                let file_name = match ttcompute_iter.next() {
                    Some(ttcompute) => format!("{idx}_{}.fastq", ttcompute.core.gene),
                    None => format!("{idx}.fastq"),
                };
                let filepath = pathbuf.join(file_name);
                myrio_core::io::write_fastq_to_file(
                    filepath,
                    cluster,
                    &myrio_core::io::CompressionMethod::None,
                )?;
            }
        }

        clusters
            .into_iter()
            .map(|myrseqs| {
                let elements = myrseqs
                    .into_iter()
                    .map(|myrseq| myrseq.compute_kmer_counts(search_k, search_t1).0)
                    .collect_vec();
                myrio_core::clustering::compute_cluster_centroid(&elements)
            })
            .collect_vec()
    };

    spinner.finish_with_message("Clustering");

    let mut csv_file_opt: Option<File> = mat
        .remove_one::<PathBuf>("output-csv")
        .map(|mut pathbuf| {
            if pathbuf.is_dir() {
                let timestamp_in_seconds =
                    SystemTime::now().duration_since(SystemTime::UNIX_EPOCH).unwrap().as_secs();
                let file_name = format!("myrio_results_{timestamp_in_seconds}.csv");
                pathbuf.set_file_name(file_name);
            } else {
                pathbuf.set_extension("csv");
            }
            OpenOptions::new().create(true).write(true).truncate(true).open(pathbuf)
        })
        .map_or(Ok(None), |v| v.map(Some))?;

    if let Some(ref mut csv_file) = csv_file_opt {
        csv_file.write_all(TaxTreeResults::CSV_HEADER.as_bytes())?;
    }

    let mut txt_file_opt: Option<File> = mat
        .remove_one::<PathBuf>("output-txt")
        .map(|mut pathbuf| {
            if pathbuf.is_dir() {
                let timestamp_in_seconds =
                    SystemTime::now().duration_since(SystemTime::UNIX_EPOCH).unwrap().as_secs();
                let file_name = format!("myrio_results_{timestamp_in_seconds}.txt");
                pathbuf.set_file_name(file_name);
            } else {
                pathbuf.set_extension("txt");
            }
            OpenOptions::new().create(true).write(true).truncate(true).open(pathbuf)
        })
        .map_or(Ok(None), |v| v.map(Some))?;

    for (ttcompute, query) in ttcompute_vec.into_iter().zip(queries) {
        let ttresults_full = TaxTreeResults::from_compute_tree(
            query,
            ttcompute,
            search_similarity,
            config.search.max_amount_of_leaves_per_branch,
            config.search.nb_best_analysis,
            config.search.lambda_leaf,
            config.search.lambda_branch,
            config.search.mu,
            config.search.gamma,
            config.search.delta,
            config.search.epsilon,
            None,
        )?;

        if let Some(ref mut csv_file) = csv_file_opt {
            csv_file.write_all(ttresults_full.generate_csv_records().as_bytes())?;
        }

        if let Some(ref mut txt_file) = txt_file_opt {
            if matches!(color_choice, ColorChoice::Auto) {
                console::set_colors_enabled(false);
                let prev = console::colors_enabled();
                txt_file.write_all((ttresults_full.core.to_string() + "\n").as_bytes())?;
                console::set_colors_enabled(prev);
            } else {
                txt_file.write_all((ttresults_full.core.to_string() + "\n").as_bytes())?;
            }
        }

        let ttresults_best = ttresults_full.cut(config.search.nb_best_display);

        println!("{}", ttresults_best.core);
    }

    if let Some(ref mut csv_file) = csv_file_opt {
        csv_file.sync_all()?;
    }

    if let Some(ref mut txt_file) = txt_file_opt {
        txt_file.sync_all()?;
    }

    Ok(())
}
