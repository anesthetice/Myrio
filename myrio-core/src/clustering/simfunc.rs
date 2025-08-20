use std::f64;

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
    NegEuclidDist,
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
            Self::NegEuclidDist => {
                todo!()
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
                let dot = a.merge_and_apply(b, |x, y| x * y).sum();
                SimScore::try_new(dot / (a.norm_l2() * b.norm_l2())).unwrap()
            }
            Self::Overlap => SimScore::try_new(a.merge_and_apply(b, |x, y| x.min(y)).sum()).unwrap(),
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
            Self::NegEuclidDist => {
                todo!()
            }
        }
    }
}

#[cfg(test)]
mod test {
    use std::f64;

    use bio_seq::prelude::*;

    use super::*;
    use crate::{assert_float_eq, data::MyrSeq};

    #[test]
    fn simfunc_test() {
        let myrseq_1 = MyrSeq::create("", None, dna!["ACTG"], &[0; 4]);
        let myrseq_2 = MyrSeq::create("", None, dna!["ACTC"], &[0; 4]);

        let a = myrseq_1.compute_dense_kmer_counts(2, f64::MIN).unwrap().0;
        let b = myrseq_2.compute_dense_kmer_counts(2, f64::MIN).unwrap().0;
        let cosine_dense = SimilarityFunction::Cosine.compute_dense(&a, &b);
        assert_float_eq!(*cosine_dense, 2.0 / 3.0);
        let overlap_dense = SimilarityFunction::OverlapNorm.compute_dense(&a, &b);
        assert_float_eq!(*overlap_dense, 0.5);

        let a = myrseq_1.compute_sparse_kmer_counts(2, f64::MIN).unwrap().0;
        let b = myrseq_2.compute_sparse_kmer_counts(2, f64::MIN).unwrap().0;
        let cosine_sparse = SimilarityFunction::Cosine.compute_sparse(&a, &b);
        assert_float_eq!(*cosine_sparse, 2.0 / 3.0);
        let overlap_sparse = SimilarityFunction::OverlapNorm.compute_sparse(&a, &b);
        assert_float_eq!(*overlap_sparse, 0.5);

        assert_float_eq!(*cosine_dense, *cosine_sparse);
        assert_float_eq!(*overlap_dense, *overlap_sparse);
    }
}
