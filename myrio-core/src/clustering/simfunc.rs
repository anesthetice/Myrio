// Imports
use itertools::Itertools;
use ndarray::Array1;

use crate::data::SFVec;
use nutype::nutype;

pub type SimScore = SimilarityScore;
pub type SimFunc = SimilarityFunction;

#[nutype(
    default = 0_f64,
    sanitize(with = |val| val.clamp(0.0, 1.0)),
    validate(greater_or_equal = 0.0, less_or_equal = 1.0, finite),
    derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deref, TryFrom, Display, Default)
)]
pub struct SimilarityScore(f64);

#[derive(Debug, Clone, Copy)]
pub enum SimilarityFunction {
    Cosine,
    Overlap,
}

impl SimilarityFunction {
    pub fn compute_dense(
        &self,
        a: &Array1<f64>,
        b: &Array1<f64>,
    ) -> SimScore {
        match &self {
            Self::Cosine => SimScore::try_new(
                a.dot(b) / (a.mapv(|val| val * val).sum().sqrt() * b.mapv(|val| val * val).sum().sqrt()),
            )
            .unwrap(),
            Self::Overlap => {
                let mut max = Array1::<f64>::zeros(a.len());
                let mut min = Array1::<f64>::zeros(a.len());

                for (idx, (&a, &b)) in a.iter().zip_eq(b.iter()).enumerate() {
                    max[idx] = a.max(b);
                    min[idx] = a.min(b);
                }

                SimilarityScore::try_new(min.sum() / max.sum()).unwrap()
            }
        }
    }

    pub fn compute_sparse(
        &self,
        a: &SFVec,
        b: &SFVec,
    ) -> SimScore {
        match &self {
            Self::Cosine => {
                let dot: f64 = {
                    if a.len() > b.len() {
                        b.iter()
                            .map(|(&key, &val)| match a.get(key) {
                                Some(&other_val) => val * other_val,
                                None => 0.0,
                            })
                            .sum()
                    } else {
                        a.iter()
                            .map(|(&key, &val)| match b.get(key) {
                                Some(&other_val) => val * other_val,
                                None => 0.0,
                            })
                            .sum()
                    }
                };
                SimScore::try_new(
                    dot / (a.values().map(|val| val * val).sum::<f64>().sqrt()
                        * b.values().map(|val| val * val).sum::<f64>().sqrt()),
                )
                .unwrap()
            }
            Self::Overlap => {
                let mut max = Array1::<f64>::zeros(a.len() + b.len());
                let mut min = Array1::<f64>::zeros(a.len() + b.len());

                for (idx, &key) in a.keys().chain(b.keys()).unique().enumerate() {
                    max[idx] = a[key].max(b[key]);
                    min[idx] = a[key].min(b[key]);
                }

                SimScore::try_new(min.sum() / max.sum()).unwrap()
            }
        }
    }
}

#[cfg(test)]
mod test {
    use std::f64;

    use bio_seq::prelude::*;

    use super::*;
    use crate::data::MyrSeq;

    #[test]
    fn simfunc_test() {
        let myrseq_1 = MyrSeq::create("", None, dna!["ACTG"], &[0; 4]);
        let myrseq_2 = MyrSeq::create("", None, dna!["ACTC"], &[0; 4]);

        let a = myrseq_1.compute_dense_kmer_counts(2, f64::MIN).unwrap().0;
        let b = myrseq_2.compute_dense_kmer_counts(2, f64::MIN).unwrap().0;
        let cosine_dense = SimilarityFunction::Cosine.compute_dense(&a, &b);
        assert!((*cosine_dense - 2.0 / 3.0).abs() < f64::EPSILON);
        let overlap_dense = SimilarityFunction::Overlap.compute_dense(&a, &b);
        assert!((*overlap_dense - 0.5).abs() < f64::EPSILON);

        let a = myrseq_1.compute_sparse_kmer_counts(2, f64::MIN).unwrap().0;
        let b = myrseq_2.compute_sparse_kmer_counts(2, f64::MIN).unwrap().0;
        let cosine_sparse = SimilarityFunction::Cosine.compute_sparse(&a, &b);
        assert!((*cosine_sparse - 2.0 / 3.0).abs() < f64::EPSILON);
        let overlap_sparse = SimilarityFunction::Overlap.compute_sparse(&a, &b);
        assert!((*overlap_sparse - 0.5).abs() < f64::EPSILON);

        assert_eq!(cosine_dense, cosine_sparse);
        assert_eq!(overlap_dense, overlap_sparse);
    }
}
