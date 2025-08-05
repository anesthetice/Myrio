// Modules
mod simfunc;

// Re-exports
pub use simfunc::SimilarityFunction;

// Imports
use crate::MyrSeq;
use itertools::Itertools;
use ndarray::Array1;

pub struct Clusterer {}

impl Clusterer {
    pub fn cluster(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1_cutoff: f64,
        t2_cutoff: f64,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        if MyrSeq::K_DENSE_VALID_RANGE.contains(&k) {
            Self::_cluster_dense(myrseqs, k, t1_cutoff, t2_cutoff, similarity_function)
        } else if MyrSeq::K_SPARSE_VALID_RANGE.contains(&k) {
            Self::_cluster_sparse(myrseqs, k, t1_cutoff, t2_cutoff, similarity_function)
        } else {
            panic!("Only k ∈ {{2, ..., 42}} is currently supported")
        }
    }

    pub fn _cluster_dense(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1_cutoff: f64,
        t2_cutoff: f64,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        struct Cluster {
            seeds: Array1<f64>,
            elements: Vec<MyrSeq>,
        }

        impl Cluster {
            fn merge(
                &mut self,
                b: Self,
            ) {
                self.seeds = (&self.seeds + &b.seeds).mapv(|val| val / 2.0);
                self.elements.reserve(b.elements.len());
                self.elements.extend(b.elements);
            }
        }

        // Step 1, collect every myrseq as a cluster
        let mut clusters: Vec<Cluster> = myrseqs
            .into_iter()
            .map(|myrseq| {
                let (map, _) = myrseq.compute_dense_kmer_counts(k, t1_cutoff).unwrap();
                Cluster { seeds: map, elements: vec![myrseq] }
            })
            .collect_vec();

        // Step 2, continuously iterate over clusters and merge as much as possible
        let mut halt = false;
        while !halt {
            halt = true;
            let mut i = 0;
            while i < clusters.len() {
                let mut j = i + 1;
                while j < clusters.len() {
                    if similarity_function.compute_dense(&clusters[i].seeds, &clusters[j].seeds) > t2_cutoff {
                        halt = false;
                        let cluster_to_be_merged = clusters.remove(j);
                        clusters[i].merge(cluster_to_be_merged);
                    }
                    j += 1;
                }
                i += 1;
            }
        }

        clusters.into_iter().map(|cluster| cluster.elements).collect_vec()
    }

    pub fn _cluster_sparse(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1_cutoff: f64,
        t2_cutoff: f64,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        todo!()
    }
}
