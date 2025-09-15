/*
pub fn clustering_assessment() -> anyhow::Result<()> {
    let testset = load_testset()?;
    let testset_ref = &testset;

    let inner = |k: usize, t1: f64, simfunc: SimFunc| {
        println!("for k={k}, t1={t1}, simfunc={simfunc:?}");
        let sum: f64 = testset_ref
            .clone()
            .into_iter()
            .map(|(myrseqs, nb_clusters, info)| {
                println!("## {info}");
                let clusters = crate::clustering::partition::Clusterer::_cluster_dense(
                    myrseqs,
                    nb_clusters,
                    k,
                    t1,
                    simfunc,
                );
                #[cfg(debug_assertions)]
                for (idx, cluster) in clusters.iter().enumerate() {
                    let mut id_count_map: HashMap<&str, usize> = HashMap::new();
                    for myrseq in cluster.iter() {
                        let count_ref = id_count_map.entry(myrseq.id.as_str()).or_default();
                        *count_ref += 1;
                    }
                    println!("cluster {idx}: {id_count_map:?}")
                }
                #[cfg(debug_assertions)]
                println!();

                compute_cluster_cost(clusters)
            })
            .sum();
        println!("→ total cost = {sum:.3}  (higher is worse)\n");
    };

    // did really well
    inner(6, 0.2, SimFunc::Cosine);
    inner(6, 0.2, SimFunc::Overlap);

    Ok(())
}
*/
