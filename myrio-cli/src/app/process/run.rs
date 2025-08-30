use indicatif::MultiProgress;
use itertools::Itertools;
use myrio_core::{
    clustering::{ClusterInitializationMethod, ClusteringParameters},
    similarity::SimFunc,
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
    const K_SEARCH: usize = 20;
    const REPR_SAMPLES: usize = 200;
    const CACHE_OPT: CacheOptions = CacheOptions::Disabled;

    const T1: f64 = 0.2;
    const T2: usize = 5;
    const SIMFUNC: SimFunc = myrio_core::similarity::cosine_similarity_already_normalized;
    const MULTITHREADING_FLAG: bool = true;

    let myrseqs: Vec<MyrSeq> = myrio_core::io::read_fastq_from_file(input_filepath)?;

    // gather compute trees
    let mut rng = rand::rngs::StdRng::try_from_os_rng()?;
    let compute_trees = tree_filepaths
        .iter()
        .map(|filepath| {
            let multi = MultiProgress::new();
            let spinner = myrio_core::utils::simple_spinner(
                Some(format!("Loading and processing '{}'", filepath.display())),
                Some(200),
                Some(&multi),
            );
            let ttcompute = TaxTreeCompute::from_store_tree(
                TaxTreeStore::decode_from_file(filepath)?,
                K_CLUSTER,
                K_SEARCH,
                REPR_SAMPLES,
                config.fasta_expansion_max_consecutive_N_before_gap,
                CACHE_OPT,
                &mut rng,
                Some(&multi),
            );
            spinner.finish();
            ttcompute
        })
        .collect::<Result<Vec<TaxTreeCompute>, myrio_core::tax::Error>>()?;

    let spinner = myrio_core::utils::simple_spinner(Some("Clustering...".to_string()), Some(200), None);

    let cim = ClusterInitializationMethod::FromCentroids(
        compute_trees.iter().map(TaxTreeCompute::get_kmer_normcounts_repr).collect_vec(),
    );
    let params = ClusteringParameters::new(K_CLUSTER, T1, T2, SIMFUNC, cim, MULTITHREADING_FLAG);

    let queries = myrio_core::clustering::cluster(myrseqs, params)
        .clusters
        .into_iter()
        .map(|myrseqs| {
            let elements = myrseqs
                .into_iter()
                .map(|myrseq| myrseq.compute_sparse_kmer_counts(K_SEARCH, T1).0)
                .collect_vec();
            myrio_core::clustering::compute_cluster_centroid_without_parallelism(&elements)
        })
        .collect_vec();

    spinner.finish_with_message("Finished clustering");

    for (ctree, query) in compute_trees.into_iter().zip_eq(queries) {
        let res_tree = TaxTreeResults::from_compute_tree(query, ctree, SIMFUNC, None)?;
        res_tree.test();
    }

    Ok(())
}
