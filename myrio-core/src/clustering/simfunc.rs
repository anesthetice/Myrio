// Imports
use itertools::Itertools;
use nutype::nutype;

use crate::data::{DFArray, SFVec};

pub type SimScore = SimilarityScore;
pub type SimFunc = SimilarityFunction;

#[nutype(
    default = 0_f64,
    validate(finite),
    derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deref, TryFrom, Display, Default)
)]
pub struct SimilarityScore(f64);

#[derive(Debug, Clone, Copy)]
pub enum SimilarityFunction {
    Cosine,
    Overlap,
    OverlapNorm,
}

impl SimilarityFunction {
    pub fn compute_dense(
        &self,
        a: &DFArray,
        b: &DFArray,
    ) -> SimScore {
        match &self {
            Self::Cosine => SimScore::try_new(
                a.dot(b) / (a.mapv(|val| val * val).sum().sqrt() * b.mapv(|val| val * val).sum().sqrt()),
            )
            .unwrap(),
            Self::Overlap => {
                let mut sum_min: f64 = 0.0;
                ndarray::Zip::from(a).and(b).for_each(|&x, &y| sum_min += x.min(y));
                SimilarityScore::try_new(sum_min).unwrap()
            }
            Self::OverlapNorm => {
                let (mut sum_min, mut sum_max) = (0.0, 0.0);

                ndarray::Zip::from(a).and(b).for_each(|&x, &y| {
                    sum_min += x.min(y);
                    sum_max += x.max(y);
                });

                SimilarityScore::try_new(sum_min / sum_max).unwrap()
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
                let dot: f64 = if a.count() > b.count() {
                    b.iter().filter_map(|(&key, &b_val)| a.get(key).map(|&a_val| a_val * b_val)).sum()
                } else {
                    a.iter().filter_map(|(&key, &a_val)| b.get(key).map(|&b_val| a_val * b_val)).sum()
                };
                SimScore::try_new(
                    dot / (a.values().map(|val| val * val).sum::<f64>().sqrt()
                        * b.values().map(|val| val * val).sum::<f64>().sqrt()),
                )
                .unwrap()
            }
            Self::Overlap => {
                let mut sum_min = 0.0;
                for key in a.merge_keys(b) {
                    sum_min += a[key].min(b[key]);
                }
                SimScore::try_new(sum_min).unwrap()
            }
            Self::OverlapNorm => {
                let (mut sum_min, mut sum_max) = (0.0, 0.0);

                for key in a.merge_keys(b) {
                    let a_val = a[key];
                    let b_val = b[key];
                    sum_min += a_val.min(b_val);
                    sum_max += a_val.max(b_val);
                }

                SimScore::try_new(sum_min / sum_max).unwrap()
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
