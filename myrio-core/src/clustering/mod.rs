// Imports
use crate::MyrSeq;

pub struct Clusterer {}

impl Clusterer {
    pub fn cluster(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1_cutoff: f64,
        t2_cutoff: f64,
        similarity_function: SimilarityFunction,
    ) {
        if MyrSeq::K_DENSE_VALID_RANGE.contains(&k) {
            Self::_cluster_dense(myrseqs, k, t1_cutoff, t2_cutoff, similarity_function);
        } else if MyrSeq::K_SPARSE_VALID_RANGE.contains(&k) {
            Self::_cluster_sparse(myrseqs, k, t1_cutoff, t2_cutoff, similarity_function);
        } else {
            ()
        }
    }

    pub fn _cluster_dense(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1_cutoff: f64,
        t2_cutoff: f64,
        similarity_function: SimilarityFunction,
    ) {
    }

    pub fn _cluster_sparse(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1_cutoff: f64,
        t2_cutoff: f64,
        similarity_function: SimilarityFunction,
    ) {
    }
}

pub enum SimilarityFunction {
    Cosine,
    Overlap,
}
