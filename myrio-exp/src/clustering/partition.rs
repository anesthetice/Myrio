// Imports
use std::{collections::HashMap, f64};

use itertools::Itertools;
use myrio_core::{
    data::{Float, MyrSeq, SFVec},
    similarity::SimScore,
};

use crate::{DFArray, compute_dense_kmer_counts, simfunc::SimilarityFunction};

/// Partition-based clustering
pub struct Clusterer;

impl Clusterer {
    pub fn _cluster_dense(
        myrseqs: Vec<MyrSeq>,
        nb_clusters: usize,
        k: usize,
        t1: Float,
        similarity_function: SimilarityFunction,
    ) -> Vec<DFArray> {
        //Vec<Vec<MyrSeq>> {
        struct Cluster<'a> {
            centroid_seeds: DFArray,
            elements: Vec<&'a DFArray>,
        }

        impl<'a> Cluster<'a> {
            fn new(seeds: &'a DFArray) -> Self {
                Self { centroid_seeds: seeds.clone(), elements: vec![seeds] }
            }

            // TODO: rename
            fn reincarnate(
                self,
                countmap_size: usize,
            ) -> Self {
                let n = self.elements.len() as f64;
                if self.elements.is_empty() {
                    return Self { centroid_seeds: self.centroid_seeds, elements: Vec::new() };
                }
                let new_centroid_seeds: DFArray = self
                    .elements
                    .into_iter()
                    .fold::<DFArray, _>(DFArray::zeros(countmap_size), |acc, x| acc + x)
                    .mapv(|v| v / n);

                Self { centroid_seeds: new_centroid_seeds, elements: Vec::new() }
            }
        }

        let countmap_size = 4_usize.pow(k as u32);

        // Step 1, compute the k-mer counts for each myrseq
        let (myrseqs, kcounts, nb_hcks): (Vec<MyrSeq>, Vec<DFArray>, Vec<usize>) = myrseqs
            .into_iter()
            .filter_map(|myrseq| {
                let (kcount, nb_hck) = compute_dense_kmer_counts(&myrseq, k, t1);
                if nb_hck != 0 { Some((myrseq, kcount, nb_hck)) } else { None }
            })
            .sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
            .multiunzip();

        // Step 2, initialize the clusters
        let mut clusters: Vec<Cluster> = vec![Cluster::new(&kcounts[0])];
        let max_nb_hck = nb_hcks[0] as Float;
        while clusters.len() < nb_clusters {
            let new_cluster_seeds = kcounts
                .iter()
                .zip_eq(nb_hcks.iter())
                .min_by_key(|&(kcount, nb_hck)| {
                    let mut score: Float = clusters
                        .iter()
                        .map(|cluster| *similarity_function.compute_dense(&cluster.centroid_seeds, kcount))
                        .sum();

                    // mean similarity score
                    score /= clusters.len() as Float;
                    // penalty (increases score) for seqs with a low number of high-confidence k-mers
                    score = 0.7 * score + 0.3 * (1.0 - (*nb_hck as Float) / max_nb_hck);
                    SimScore::try_new(score).unwrap()
                })
                .unwrap()
                .0;
            clusters.push(Cluster::new(new_cluster_seeds));
        }

        // Step 3
        for _ in 0..5 {
            for kcount in kcounts.iter() {
                clusters
                    .iter_mut()
                    .max_by_key(|cluster| similarity_function.compute_dense(&cluster.centroid_seeds, kcount))
                    .unwrap()
                    .elements
                    .push(kcount);
            }

            clusters = clusters.into_iter().map(|cluster| cluster.reincarnate(countmap_size)).collect_vec();
        }

        clusters.into_iter().map(|cl| cl.centroid_seeds).collect()

        /*
        // Step 4
        let mut output: Vec<Vec<MyrSeq>> = vec![Vec::new(); nb_clusters];
        for (myrseq, kcount) in myrseqs.into_iter().zip_eq(kcounts.iter()) {
            let idx = clusters
                .iter()
                .enumerate()
                .max_by_key(|(idx, cluster)| {
                    similarity_function.compute_dense(&cluster.centroid_seeds, kcount)
                })
                .unwrap()
                .0;
            output[idx].push(myrseq);
        }

        output
        */
    }

    pub fn _cluster_sparse(
        myrseqs: Vec<MyrSeq>,
        nb_clusters: usize,
        k: usize,
        t1: Float,
        similarity_function: SimilarityFunction,
    ) -> Vec<SFVec> {
        struct Cluster<'a> {
            centroid_seeds: SFVec,
            elements: Vec<&'a SFVec>,
        }

        impl<'a> Cluster<'a> {
            fn new(seeds: &'a SFVec) -> Self {
                Self { centroid_seeds: seeds.clone(), elements: vec![seeds] }
            }

            // TODO: rename
            fn reincarnate(
                self,
                countmap_size: usize,
            ) -> Self {
                if self.elements.is_empty() {
                    return Self { centroid_seeds: self.centroid_seeds, elements: Vec::new() };
                }

                let n = self.elements.len() as Float;
                let new_centroid_seeds: SFVec =
                    self.elements.into_iter().fold::<SFVec, _>(SFVec::new(countmap_size), |acc, x| {
                        acc.merge_and_apply(x, |a, b| a + b)
                    }) / n;

                Self { centroid_seeds: new_centroid_seeds, elements: Vec::new() }
            }
        }

        let countmap_size = 4_usize.pow(k as u32);

        // Step 1, compute the k-mer counts for each myrseq
        let (myrseqs, kcounts, nb_hcks): (Vec<MyrSeq>, Vec<SFVec>, Vec<usize>) = myrseqs
            .into_iter()
            .filter_map(|myrseq| {
                let (mut kcount, nb_hck) = myrseq.compute_kmer_counts(k, t1);
                if nb_hck != 0 { Some((myrseq, kcount, nb_hck)) } else { None }
            })
            .sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
            .multiunzip();

        // Step 2, initialize the clusters
        let mut clusters: Vec<Cluster> = vec![Cluster::new(&kcounts[0])];
        let max_nb_hck = nb_hcks[0] as Float;
        while clusters.len() < nb_clusters {
            let new_cluster_seeds = kcounts
                .iter()
                .zip_eq(nb_hcks.iter())
                .min_by_key(|&(kcount, nb_hck)| {
                    let mut score: Float = clusters
                        .iter()
                        .map(|cluster| *similarity_function.compute_sparse(&cluster.centroid_seeds, kcount))
                        .sum();

                    // mean similarity score
                    score /= clusters.len() as Float;
                    // penalty (increases score) for seqs with a low number of high-confidence k-mers
                    score = 0.7 * score + 0.3 * (1.0 - (*nb_hck as Float) / max_nb_hck);
                    SimScore::try_new(score).unwrap()
                })
                .unwrap()
                .0;
            clusters.push(Cluster::new(new_cluster_seeds));
        }

        // Step 3
        for _ in 0..3 {
            for kcount in kcounts.iter() {
                clusters
                    .iter_mut()
                    .max_by_key(|cluster| similarity_function.compute_sparse(&cluster.centroid_seeds, kcount))
                    .unwrap()
                    .elements
                    .push(kcount);
            }

            clusters = clusters.into_iter().map(|cluster| cluster.reincarnate(countmap_size)).collect_vec();
        }
        clusters.into_iter().map(|cl| cl.centroid_seeds).collect()
        /*
        // Step 4
        let mut output: Vec<Vec<MyrSeq>> = vec![Vec::new(); nb_clusters];
        for (myrseq, kcount) in myrseqs.into_iter().zip_eq(kcounts.iter()) {
            let idx = clusters
                .iter()
                .enumerate()
                .max_by_key(|(idx, cluster)| {
                    similarity_function.compute_sparse(&cluster.centroid_seeds, kcount)
                })
                .unwrap()
                .0;
            output[idx].push(myrseq);
        }

        output
        */
    }
}

pub fn compute_cluster_cost(clusters: Vec<Vec<MyrSeq>>) -> f64 {
    let mut cost: f64 = 0.0;

    for (idx, cluster) in clusters.into_iter().enumerate() {
        let mut id_count_map: HashMap<&str, usize> = HashMap::new();
        for myrseq in cluster.iter() {
            let count_ref = id_count_map.entry(myrseq.id.as_str()).or_default();
            *count_ref += 1;
        }
        if !id_count_map.is_empty() {
            let counts = id_count_map.values().copied().collect_vec();
            let max_count = *counts.iter().max().unwrap();
            let diff = (counts.iter().sum::<usize>() - max_count) as f64;
            cost += diff
        }
    }
    cost
}
