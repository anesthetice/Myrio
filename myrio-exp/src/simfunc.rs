// Imports
use std::f64;

use myrio_core::{
    data::{DFArray, SFVec},
    similarity::SimScore,
};

pub type SimFunc = SimilarityFunction;

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
                SimScore::try_new(sum_min).unwrap()
            }
            Self::OverlapNorm => {
                let (mut sum_min, mut sum_max) = (0.0, 0.0);

                ndarray::Zip::from(a).and(b).for_each(|&x, &y| {
                    sum_min += x.min(y);
                    sum_max += x.max(y);
                });

                SimScore::try_new(sum_min / sum_max).unwrap()
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
                let dot = a.merge_with_and_apply(b, |x, y| x * y).sum();
                SimScore::try_new(dot / (a.norm_l2() * b.norm_l2())).unwrap()
            }
            Self::Overlap => SimScore::try_new(a.merge_with_and_apply(b, |x, y| x.min(y)).sum()).unwrap(),
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
