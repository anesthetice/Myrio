// Imports
use itertools::Itertools;
use myrio_core::{
    constants::Q_TO_BP_CALL_CORRECT_PROB_MAP,
    data::{Float, MyrSeq},
};

use crate::{DFArray, compute_dense_kmer_counts, simfunc::SimilarityFunction};

pub struct Clusterer;

impl Clusterer {
    pub fn cluster(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1_cutoff: Float,
        t2_cutoff: Float,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        if crate::K_DENSE_VALID_RANGE.contains(&k) {
            Self::_cluster_dense(myrseqs, k, t1_cutoff, t2_cutoff, similarity_function)
        } else if MyrSeq::K_SPARSE_VALID_RANGE.contains(&k) {
            Self::_cluster_sparse(myrseqs, k, t1_cutoff, t2_cutoff, similarity_function)
        } else {
            panic!("Only k ∈ {{2, ..., 32}} is currently supported")
        }
    }

    pub fn cluster_test(
        myrseqs: Vec<MyrSeq>,
        k: usize,
    ) -> Vec<Vec<MyrSeq>> {
        let probs = myrseqs
            .iter()
            .map(|myrseq| {
                myrseq.quality.iter().map(|q| Q_TO_BP_CALL_CORRECT_PROB_MAP[*q as usize]).collect_vec()
            })
            .concat();

        let prob_mean = probs.iter().sum::<Float>() / probs.len() as Float;
        let prob_std = (probs.iter().map(|&p| (p - prob_mean) * (p - prob_mean)).sum::<Float>()
            / (probs.len() - 1) as Float)
            .sqrt();

        println!("mean: {prob_mean:.2}\nstd: {prob_std:.2}");
        let t1_cutoff = prob_mean - 0.4 * prob_std;
        let t2_cutoff = 0.5;
        println!("T1={t1_cutoff:.2}\nT2={t2_cutoff:.2}");

        Self::cluster(myrseqs, k, t1_cutoff, t2_cutoff, SimilarityFunction::Cosine)
    }

    pub fn _cluster_dense(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1_cutoff: Float,
        t2_cutoff: Float,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        struct Cluster {
            seeds: DFArray,
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
                let (map, _nb) = compute_dense_kmer_counts(&myrseq, k, t1_cutoff);
                //println!("number of k-mers kept: {_nb} / {}", myrseq.sequence.len() - k + 1);
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
                    if *similarity_function.compute_dense(&clusters[i].seeds, &clusters[j].seeds) > t2_cutoff
                    {
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
        t1_cutoff: Float,
        t2_cutoff: Float,
        similarity_function: SimilarityFunction,
    ) -> Vec<Vec<MyrSeq>> {
        unimplemented!()
    }
}

/*
/// Based off the clustering method used by `isONclust3`; notably, we use k-mers instead of minimizers as our expected amplicon length isn't very high (<2000 bp).
pub fn isONcluster(
    myrseqs: Vec<MyrSeq>,
    k: usize,
    t1_cutoff: f64,
    t2_cutoff: f64,
    similarity_function: DenseSimFunc,
) -> anyhow::Result<Vec<Vec<MyrSeq>>> {
    if !MyrSeq::K_SPARSE_VALID_RANGE.contains(&k) {
        bail!(MyrSeq::K_SPARSE_VALID_RANGE_ERROR_MSG);
    }
    /*
    // Step 1
    let mut myrseqs_extra = myrseqs
        .into_iter()
        .map(|myrseq| {
            let (map, nb_hck) = myrseq.compute_dense_kmer_counts(k, t1_cutoff).unwrap();
            (myrseq, map, nb_hck)
        })
        .sorted_by(|(.., a), (.., b)| a.cmp(b)) // ascending order
        .collect_vec();

    // Step 2
    struct Cluster {
        reference_seeds: Vec<f64>,
        elements: Vec<MyrSeq>,
    }

    let mut clusters: Vec<Cluster> = Vec::new();
    {
        let (myrseq, ref_seeds, _) = myrseqs_extra.pop().unwrap();
        clusters.push(Cluster { reference_seeds: ref_seeds, elements: vec![myrseq] });
    }

    while let Some((myrseq, seeds, _)) = myrseqs_extra.pop() {
        match clusters
            .iter_mut()
            .find(|cluster| similarity_function(&seeds, &cluster.reference_seeds) > t2_cutoff)
        {
            Some(cluster) => cluster.elements.push(myrseq),
            None => clusters.push(Cluster { reference_seeds: seeds, elements: vec![myrseq] }),
        };
    }

    Ok(clusters.into_iter().map(|cluster| cluster.elements).collect_vec())
    */
    unimplemented!()
}
*/
