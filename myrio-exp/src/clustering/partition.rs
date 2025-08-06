use std::f64;

use itertools::Itertools;
// Imports
use myrio_core::{
    clustering::{SimScore, SimilarityFunction},
    data::MyrSeq,
};
use ndarray::Array1;

/// Partition-based clustering
pub struct Clusterer;

impl Clusterer {
    pub fn cluster(
        myrseqs: Vec<MyrSeq>,
        nb_clusters: usize,
        k: usize,
        t1: f64,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        if MyrSeq::K_DENSE_VALID_RANGE.contains(&k) {
            Self::_cluster_dense(myrseqs, nb_clusters, k, t1, similarity_function)
        } else if MyrSeq::K_SPARSE_VALID_RANGE.contains(&k) {
            Self::_cluster_sparse(myrseqs, k, nb_clusters, similarity_function)
        } else {
            panic!("Only k ∈ {{2, ..., 42}} is currently supported")
        }
    }

    pub fn _cluster_dense(
        myrseqs: Vec<MyrSeq>,
        nb_clusters: usize,
        k: usize,
        t1: f64,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        struct Cluster<'a> {
            centroid_seeds: Array1<f64>,
            elements: Vec<&'a Array1<f64>>,
        }

        impl<'a> Cluster<'a> {
            fn new(seeds: &'a Array1<f64>) -> Self {
                Self { centroid_seeds: seeds.clone(), elements: vec![seeds] }
            }
            fn reincarnate(
                self,
                countmap_size: usize,
            ) -> Self {
                let n = self.elements.len() as f64;
                let new_centroid_seeds: Array1<f64> = self
                    .elements
                    .into_iter()
                    .fold::<Array1<f64>, _>(Array1::zeros(countmap_size), |acc, x| acc + x)
                    .mapv(|v| v / n);

                Self { centroid_seeds: new_centroid_seeds, elements: Vec::new() }
            }
        }

        let countmap_size = 4_usize.pow(k as u32);

        // Step 1, compute the k-mer counts for each myrseq
        let (myrseqs, kcounts, nb_hcks): (Vec<MyrSeq>, Vec<Array1<f64>>, Vec<usize>) = myrseqs
            .into_iter()
            .map(|myrseq| {
                let (kcount, nb_hck) = myrseq.compute_dense_kmer_counts(k, t1).unwrap();
                (myrseq, kcount, nb_hck)
            })
            .sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
            .multiunzip();

        // Step 2, initialize the clusters
        let mut clusters: Vec<Cluster> = vec![Cluster::new(&kcounts[0])];
        let max_nb_hck = nb_hcks[0] as f64;
        while clusters.len() < nb_clusters {
            let new_cluster_seeds = kcounts
                .iter()
                .zip_eq(nb_hcks.iter())
                .min_by_key(|&(kcount, nb_hck)| {
                    let mut score: f64 = clusters
                        .iter()
                        .map(|cluster| *similarity_function.compute_dense(&cluster.centroid_seeds, kcount))
                        .sum();

                    // mean similarity score
                    score /= clusters.len() as f64;
                    // penalty (increases score) for seqs with a low number of high-confidence k-mers
                    score = 0.7 * score + 0.3 * (1.0 - (*nb_hck as f64) / max_nb_hck);
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
    }

    pub fn _cluster_sparse(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        nb_clusters: usize,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        todo!()
    }
}
