use std::{
    fs::{File, OpenOptions},
    io::Write,
    time::SystemTime,
};

use myrio_core::{
    clustering::ClusteringParameters,
    data::Float,
    similarity::{SimFunc, SimScore, Similarity},
    tax::{
        clade::Rank,
        compute::{CacheOptions, TaxTreeCompute},
        results::TaxTreeResults,
    },
};

use super::*;

pub fn subcommand() -> Command {
    Command::new("run")
        .arg(
            Arg::new("input")
                .help("The `.fastq` file to use as input")
                .required(true)
                //.num_args(1..101)
                .short('i')
                .visible_short_alias('q')
                .long("input")
                .visible_alias("query")
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("trees")
                .help("The one or more `.myrtree` reference databases")
                .long_help("The one or more `.myrtree` reference databases, also accepts directories")
                .required(true)
                .num_args(1..101)
                .short('t')
                .visible_short_alias('r')
                .long("trees")
                .visible_aliases(["refs", "references", "db"])
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("k-search")
                .help("The length of each k-mer (i.e., `k` itself) used for sequence comparison")
                .required(false)
                .short('k')
                .long("k-search")
                .value_parser(vparser!(usize))
                .action(ArgAction::Set)
        )
        .arg(
            Arg::new("save-clusters")
                .help("Save clusters to the specified path")
                .required(false)
                .short('s')
                .long("save-clusters")
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set)
        )
        .arg(
            Arg::new("output-csv")
                .help("Write results to a `.csv` file")
                .long_help("Write results to a `.csv` file (e.g., `--csv .` will write to a timestamped file in the current directory)")
                .required(false)
                .long("output-csv")
                .visible_aliases(["csv", "csv-output"])
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set)
        )
        .arg(
            Arg::new("output-txt")
                .help("Write results to a `.txt` file")
                .long_help("Write results to a `.txt` file (e.g., `--txt .` will write to a timestamped file in the current directory)")
                .required(false)
                .long("output-txt")
                .visible_aliases(["txt", "txt-output"])
                .value_parser(vparser!(PathBuf))
                .action(ArgAction::Set)
        )
        .arg(
            Arg::new("cache-counts")
                .help("Flag that decides if newly-computed kmer counts are then cached")
                .required(false)
                .long("cache-counts")
                .action(ArgAction::SetTrue)
        )
        .arg(
            Arg::new("nb-clusters")
                .help("The number of clusters to expect")
                .long_help("The number of clusters to expect, defaults to the number of `.myrtree` files found")
                .required(false)
                .short('n')
                .long("nb-clusters")
                .value_parser(vparser!(usize))
                .action(ArgAction::Set)
        )
}

pub fn process(
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
        let (mean_qual_score, mean_seq_len) = {
            let mut n: usize = 0;
            let mut sum: usize = 0;
            output.iter().for_each(|myrseq| {
                n += myrseq.quality.len();
                myrseq.quality.iter().for_each(|q| sum += usize::from(*q));
            });
            (sum as Float / n as Float, n as Float / output.len() as Float)
        };
        let initial_count = output.len();
        output = MyrSeq::pre_process(
            output,
            config.fastq_min_length,
            config.fastq_min_mean_qual,
            config.fastq_max_qual,
        );
        let final_count = output.len();
        printdoc! {"
            Pre-processed {initial_count} FASTQ records → ignored {} records
            {}
            {}

            ",
            initial_count - final_count,
            style(format!("├─── μ_quality_score    =  {mean_qual_score:.2}")).dim(),
            style(format!("└─── μ_sequence_length  =  {mean_seq_len:.2}")).dim(),
        };
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

            /*
            // For debugging/testing purposes
            if let Ok(ref mut ttc) = ttcompute {
                ttc.core.excise("Ficus", Rank::Genus);
            }
            */
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

    // Not the clearest name, represents for each gene, a vector of genuses (name, pool_score, conf_opt) sorted in descending order by pool_score.
    let mut ranked_genus_vec_per_gene = Vec::new();

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

        let genus_name_pool_score_conf_opt_vec = ttresults_full
            .core
            .gather_branches_at_rank(Rank::Genus)
            .into_iter()
            .map(|branch| (branch.name.clone(), branch.extra.pool_score, branch.extra.confidence))
            .sorted_unstable_by_key(|(_, score, _)| SimScore::try_new(-*score).unwrap())
            .collect_vec();

        ranked_genus_vec_per_gene.push(genus_name_pool_score_conf_opt_vec);

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

    #[inline]
    fn softmax_pooling(
        scores: Vec<Float>,
        lambda: Float,
    ) -> Float {
        let mut numerator: Float = 0.0;
        let mut denominator: Float = 0.0;

        for s in scores.into_iter() {
            numerator += s * (lambda * s).exp();
            denominator += (lambda * s).exp();
        }

        numerator / denominator
    }

    let mut top_genus_mismatch_score: Float = 0.0;
    let mut top_genus_conf_vec: Vec<Float> = Vec::new();
    for i in 0..ranked_genus_vec_per_gene.len() {
        let top_genus_name = &ranked_genus_vec_per_gene[i][0].0;
        let top_genus_conf = ranked_genus_vec_per_gene[i][0].2.unwrap().unwrap_or_else(|| {
            eprintln!("Warning, got single genus unexpectedly, assuming confidence score of 0.5");
            0.5
        });
        top_genus_conf_vec.push(top_genus_conf);
        static TAKE: usize = 8;
        for other in ranked_genus_vec_per_gene.iter().skip(i + 1) {
            let penalty = other
                .iter()
                .take(TAKE)
                .position(|(genus_name, ..)| genus_name == top_genus_name)
                .unwrap_or(TAKE) as Float
                / TAKE as Float;
            top_genus_mismatch_score += penalty * top_genus_conf.powf(0.7).max(0.12);
        }
    }

    top_genus_mismatch_score /= ranked_genus_vec_per_gene.len() as Float;

    let lambda = 0.5 + 2.0 * top_genus_conf_vec.iter().sum::<Float>();
    let top_genus_conf_score = 1.0 - softmax_pooling(top_genus_conf_vec, lambda);

    {
        let x = top_genus_mismatch_score;
        let y = top_genus_conf_score.powf(0.8);

        // Ranges from 0.0 to 1.0, a higher value indicates the sample is more likely not in the database
        // https://www.desmos.com/3d/htoedxqubn
        let unknown_score: Float = (1.0 + Float::exp(-8.0 * ((x + y) / 2.0 - 0.45))).powi(-1);

        let in_database_confidence_score = 1.0 - unknown_score;

        let colorizer = myrio_core::utils::colorizer_from_confidence_score(in_database_confidence_score);

        printdoc! {"
            In-database-confidence-score = {}
            {}
        ",
        colorizer(style(format!("{in_database_confidence_score:.3}"))),
        style("Note: experimental confidence metric, not a probability value; results are currently limited to a genus-level evaluation and may be unreliable.").dim(),
        };
    }

    Ok(())
}
