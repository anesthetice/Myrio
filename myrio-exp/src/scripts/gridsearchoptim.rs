use std::io::Write;

use itertools::Itertools;
use myrio_core::{
    clustering::ClusteringParameters,
    data::{Float, MyrSeq, SFVec},
    similarity::{SimFunc, SimScore, Similarity},
    tax::{
        clade::Rank,
        compute::{CacheOptions, TaxTreeCompute},
        results::TaxTreeResults,
        store::TaxTreeStore,
    },
};
use rand::SeedableRng;

#[rustfmt::skip]
const INPUTS: [(&str, &str); 8] = [
    ("Berberis Julianae", "./ignore/queries/Berberis_Julianae_Campus_matk_rbcL_psbA-trnH_ITS_barcode6.fastq"),
    ("Corylus Avellana", "./ignore/queries/Corylus_avellana_matk_rbcL_psbA-trnH_ITS_barcode8.fastq"),
    ("Cymbopogon Citratus", "./ignore/queries/Cymbopogon_Citratus_Qiagen_matk_rbcL_psbA-trnH_ITS_barcode10.fastq"),
    ("Echeveria Agavoides", "./ignore/queries/Echeveria_Agavoides_Fulvia_matk_rbcL_psbA-trnH_ITS_barcode4.fastq"),
    ("Echeveria Agavoides", "./ignore/queries/Echeveria_Agavoides_MN_rbcL_psbA-trnH_ITS_barcode1.fastq"),
    ("Ficus Benjamina", "./ignore/queries/Ficus_Benjamina_MN_matk_rbcL_psbA-trnH_ITS_barcode4.fastq"),
    ("Ficus Benjamina", "./ignore/queries/Ficus_Benjamina_Qiagen_matk_rbcL_psbA-trnH_ITS_barcode6.fastq"),
    ("Hedera Helix", "./ignore/queries/Hedera_Helix_Fulvia_matk_rbcL_psbA-trnH_ITS_barcode3.fastq"),
];

const TREE_FILEPATHS: [&str; 4] = [
    "./ignore/myrio-db/BOLD_Plantae_20250831_ITS.myrtree",
    "./ignore/myrio-db/BOLD_Plantae_20250831_matK.myrtree",
    "./ignore/myrio-db/BOLD_Plantae_20250831_rbcL.myrtree",
    "./ignore/myrio-db/BOLD_Plantae_20250831_trnH-psbA.myrtree",
];

const NB_CLUSTERS: usize = 4;
const NB_ITERS_MAX: usize = 10;
const MULTITHREADING_FLAG: bool = true;
const REPR_SAMPLES: usize = 200;
const CACHE_OPT: CacheOptions = CacheOptions::Disabled;

pub fn grid_search_optimization() -> anyhow::Result<()> {
    let expected_species_and_myrseqs_vec = INPUTS
        .into_iter()
        .map(|(expected_species, filepath)| {
            (expected_species, myrio_core::io::read_fastq_from_file(filepath).unwrap())
        })
        .collect_vec();

    let cluster_k_values: Vec<usize> = vec![5, 6, 7];
    let cluster_t1_values: Vec<f64> = vec![0.2];
    let cluster_t2_values: Vec<usize> = vec![20];
    let cluster_similarity_values: Vec<Similarity> =
        vec![Similarity::Cosine, Similarity::JacardTanimoto, Similarity::Overlap];
    let cluster_eta_improvement_values: Vec<Float> = vec![1E-4];
    let cluster_ssdcf_values: Vec<Float> = vec![0.4, 0.8, 1.2, 1.6, 2.0, 100.0];

    let search_k_values: Vec<usize> = vec![16, 18, 20, 22, 24, 28];
    let search_similarity_values: Vec<Similarity> =
        vec![Similarity::Cosine, Similarity::JacardTanimoto, Similarity::Overlap];

    let mut file = std::fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .open("grid_search_optimization_results.md")?;

    itertools::iproduct!(
        cluster_k_values,
        cluster_t1_values,
        cluster_t2_values,
        cluster_similarity_values,
        cluster_eta_improvement_values,
        cluster_ssdcf_values,
        search_k_values,
        search_similarity_values
    )
    .enumerate()
    .for_each(
        |(
            id,
            (
                cluster_k,
                cluster_t1,
                cluster_t2,
                cluster_similarity,
                cluster_eta_improvement_value,
                cluster_ssdcf,
                search_k,
                search_similarity,
            ),
        )| {
            let print = format!("## Entry {id}\nclustering parameters: (k={cluster_k}, t1={cluster_t1}, t2={cluster_t2}, sim={cluster_similarity:?}, η={cluster_eta_improvement_value:.4E}, ssdcf={cluster_ssdcf:.3})\nsearch parameters: (k={search_k}, sim={search_similarity:?})");
            let _ = file.write_all(print.as_bytes()).inspect_err(|e| {eprintln!("{e}")});
            println!("{print}");

            let results = single(cluster_k, cluster_t1, cluster_t2, cluster_similarity, cluster_eta_improvement_value, cluster_ssdcf, search_k, search_similarity);

            match results {
                Ok(output) => {
                    println!("{output}");
                    let _ = file.write_all(output.as_bytes()).inspect_err(|e| {eprintln!("{e}")});
                },
                Err(e) => {
                    eprintln!("Error occured: {e}\n\n");
                    let _ = writeln!(&mut file, "Error occured: {e}\n\n").inspect_err(|e| {eprintln!("{e}")});
                }
            }


        }
    );

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn single(
    //
    cluster_k: usize,
    cluster_t1: f64,
    cluster_t2: usize,
    cluster_similarity: Similarity,
    cluster_eta_improvement_value: Float,
    cluster_ssdcf: Float,
    //
    search_k: usize,
    search_similarity: Similarity,
) -> anyhow::Result<String> {
    let cluster_simfunc = cluster_similarity.to_simfunc(true);

    let cluster_params = ClusteringParameters {
        k: cluster_k,
        t1: cluster_t1,
        t2: cluster_t2,
        similarity: cluster_similarity,
        intial_centroids: Vec::new(),
        expected_nb_of_clusters: 4,
        eta_improvement: cluster_eta_improvement_value,
        nb_iters_max: NB_ITERS_MAX,
        silhouette_std_deviation_cutoff_factor: cluster_ssdcf,
    };

    // gather compute trees
    let mut rng = rand::rngs::StdRng::try_from_os_rng()?;

    let ttcompute_vec = TREE_FILEPATHS
        .iter()
        .map(|filepath| {
            TaxTreeCompute::from_store_tree(
                TaxTreeStore::decode_from_file(filepath)?,
                cluster_k,
                search_k,
                REPR_SAMPLES,
                3,
                CACHE_OPT,
                &mut rng,
                None,
            )
        })
        .collect::<Result<Vec<TaxTreeCompute>, myrio_core::tax::Error>>()?;

    let mut total_score: f64 = 0.0;
    let mut output: String = String::new();

    for (expected_species, filepath) in INPUTS.into_iter() {
        let expected_genus = expected_species.split_once(" ").unwrap().0;
        let myrseqs = myrio_core::io::read_fastq_from_file(filepath)?;

        output.push_str(expected_species);
        output.push_str(": ");

        #[rustfmt::skip]
        let (mut search_queries, mut match_queries): (Vec<SFVec>, Vec<SFVec>) =
            myrio_core::clustering::cluster(myrseqs, cluster_params.clone())
                .clusters
                .into_iter()
                .map(|myrseqs| {
                    let search_elements = myrseqs
                        .iter()
                        .map(|myrseq| myrseq.compute_sparse_kmer_counts(search_k, cluster_t1).0)
                        .collect_vec();
                    let search_queries =
                        myrio_core::clustering::compute_cluster_centroid(&search_elements);

                    let match_elements = myrseqs
                        .iter()
                        .map(|myrseq| myrseq.compute_sparse_kmer_counts(cluster_k, cluster_t1).0)
                        .collect_vec();
                    let match_queries =
                        myrio_core::clustering::compute_cluster_centroid(&match_elements);

                    (search_queries, match_queries)
                })
                .multiunzip();

        for ttcompute in ttcompute_vec.iter() {
            let query_idx = match_queries
                .iter()
                .enumerate()
                .max_by_key(|(_, q)| cluster_simfunc(ttcompute.get_kmer_normalized_counts_fingerprint(), q))
                .unwrap()
                .0;

            let query = search_queries.remove(query_idx);
            match_queries.remove(query_idx);

            let ttresults_best =
                TaxTreeResults::from_compute_tree(query, ttcompute.clone(), search_similarity, None)?.cut(5);

            let best_genus = &*ttresults_best
                .core
                .gather_branches_at_rank(Rank::Genus)
                .iter()
                .max_by_key(|branch| SimScore::try_from(branch.extra.mean).unwrap_or_default())
                .unwrap()
                .name;

            let best_species_vec = ttresults_best
                .core
                .gather_leaves()
                .into_iter()
                .sorted_by_key(|leaf| ttresults_best.core.payloads[leaf.payload_id])
                .map(|leaf| &*leaf.name)
                .collect_vec();

            println!("{}{best_genus}\n{best_species_vec:?}\n", ttresults_best.core);

            let mut score: f64 = 0.0;

            output += "(";
            output += &ttresults_best.core.gene;
            output += " ";

            if best_genus == expected_genus {
                score += 1.0;
                output.push_str("✅g ");
            } else {
                output.push_str("❌g ");
            }

            if best_species_vec[0] == expected_species {
                score += 2.0;
                output.push_str("👑s) ");
            } else if best_species_vec[0..3].contains(&expected_species) {
                score += 1.0;
                output.push_str("3️⃣s) ");
            } else if best_species_vec.contains(&expected_genus) {
                score += 0.5;
                output.push_str("5️⃣s) ");
            } else {
                output.push_str("❌s) ");
            }

            total_score += score;
        }
        output.push('\n');
        println!("{output}");
    }

    output += &format!("=> score = {total_score}\n\n");

    Ok(output)
}
