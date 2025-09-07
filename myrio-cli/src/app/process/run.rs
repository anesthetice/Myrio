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

    const K_CLUSTER: usize = 6;
    const K_SEARCH: usize = 18;
    const REPR_SAMPLES: usize = 500;
    const CACHE_OPT: CacheOptions = CacheOptions::Disabled;

    const T1: Float = 0.2;
    const T2: usize = 5;
    const SIMILARITY: Similarity = Similarity::JacardTanimoto;

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
            let ttcompute = TaxTreeCompute::from_store_tree(
                TaxTreeStore::load_from_file(filepath)?,
                K_CLUSTER,
                K_SEARCH,
                REPR_SAMPLES,
                config.nb_bootstrap_resamples,
                CACHE_OPT,
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

    let intial_centroids = if mat.get_flag("no-initial-centroids") {
        Vec::with_capacity(0)
    } else {
        ttcompute_vec
            .iter()
            .map(TaxTreeCompute::get_kmer_normalized_counts_fingerprint)
            .cloned()
            .collect_vec()
    };

    let simfunc: SimFunc = SIMILARITY.to_simfunc(true);

    let params = ClusteringParameters {
        k: K_CLUSTER,
        t1: T1,
        t2: T2,
        similarity: SIMILARITY,
        intial_centroids,
        expected_nb_of_clusters: nb_clusters,
        eta_improvement: 1E-4,
        nb_iters_max: 20,
        silhouette_std_deviation_cutoff_factor: 1.6,
    };

    // The nice "default way" where we can use cluster centroids
    if !mat.get_flag("no-initial-centroids") {
        let queries = myrio_core::clustering::cluster(myrseqs, params)
            .clusters
            .into_iter()
            .map(|myrseqs| {
                let elements = myrseqs
                    .into_iter()
                    .map(|myrseq| myrseq.compute_kmer_counts(K_SEARCH, T1).0)
                    .collect_vec();
                myrio_core::clustering::compute_cluster_centroid(&elements)
            })
            .collect_vec();

        spinner.finish_with_message("Finished clustering");

        for (ttcompute, query) in ttcompute_vec.into_iter().zip(queries) {
            let ttresults_full =
                TaxTreeResults::from_compute_tree(query, ttcompute.clone(), SIMILARITY, None)?;

            /*
            use std::io::Write;
            let mut file =
                std::fs::OpenOptions::new().write(true).create(true).truncate(true).open("./output.txt")?;
            writeln!(file, "{ttresults_full}")?;
            file.sync_all()?;
            drop(file);
            break;
            */

            let ttresults_best = ttresults_full.cut(10);

            println!("{}", ttresults_best.core);
        }
    } else {
        #[rustfmt::skip]
        let (mut search_queries, mut match_queries): (Vec<SFVec>, Vec<SFVec>) =
            myrio_core::clustering::cluster(myrseqs, params.clone())
                .clusters
                .into_iter()
                .map(|myrseqs| {
                    let search_elements = myrseqs
                        .iter()
                        .map(|myrseq| myrseq.compute_kmer_counts(K_SEARCH, T1).0)
                        .collect_vec();
                    let search_queries =
                        myrio_core::clustering::compute_cluster_centroid(&search_elements);

                    let match_elements = myrseqs
                        .iter()
                        .map(|myrseq| myrseq.compute_kmer_counts(K_CLUSTER, T1).0)
                        .collect_vec();
                    let match_queries =
                        myrio_core::clustering::compute_cluster_centroid(&match_elements);

                    (search_queries, match_queries)
                })
                .multiunzip();

        spinner.finish_with_message("Finished clustering");

        for ttcompute in ttcompute_vec.iter() {
            let query_idx = match_queries
                .iter()
                .enumerate()
                .max_by_key(|(_, q)| {
                    let score = simfunc(ttcompute.get_kmer_normalized_counts_fingerprint(), q);
                    println!("{score}");
                    score
                })
                .unwrap()
                .0;

            let query = search_queries.remove(query_idx);
            match_queries.remove(query_idx);

            let ttresults_best =
                TaxTreeResults::from_compute_tree(query, ttcompute.clone(), SIMILARITY, None)?.cut(5);

            println!("{}", ttresults_best.core);
        }
    }

    Ok(())
}
