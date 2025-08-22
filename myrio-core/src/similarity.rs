// Imports
use std::f64;

use itertools::Itertools;
use nutype::nutype;

use crate::data::SFVec;

pub type SimFunc = fn(&SFVec, &SFVec) -> SimScore;

#[nutype(
    default = 0_f64,
    validate(finite),
    derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deref, TryFrom, Display, Default)
)]
pub struct SimScore(f64);

pub fn cosine_similarity(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.merge_with_and_apply(b, |x, y| x * y).sum();
    (dot / (a.norm_l2() * b.norm_l2())).try_into().unwrap()
}

pub fn cosine_similarity_already_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    a.merge_with_and_apply(b, |x, y| x * y).sum().try_into().unwrap()
}

pub fn overlap_similarity(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let (_, min_max, dim, (min_s, max_s)) = a.merge_and_apply(b, |x, y| (x.min(y), x.max(y)));

    let nb_sparse = dim - min_max.len();
    let sum_min = min_max.iter().map(|(min, _)| min).sum1().unwrap_or(0.0) + min_s * nb_sparse as f64;
    let sum_max = min_max.iter().map(|(_, max)| max).sum1().unwrap_or(0.0) + max_s * nb_sparse as f64;

    SimScore::try_new(sum_min / sum_max).unwrap()
}

pub fn overlap_similarity_already_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    a.merge_with_and_apply(b, |x, y| x.min(y)).sum().try_into().unwrap()
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

        let mut a = myrseq_1.compute_sparse_kmer_counts(2, f64::MIN).unwrap().0;
        let mut b = myrseq_2.compute_sparse_kmer_counts(2, f64::MIN).unwrap().0;
        let cosine = cosine_similarity(&a, &b);
        assert_float_eq!(*cosine, 2.0 / 3.0);
        let overlap = overlap_similarity(&a, &b);
        assert_float_eq!(*overlap, 0.5);

        a.normalize_l2();
        b.normalize_l2();

        assert_float_eq!(*cosine_similarity_already_normalized(&a, &b), *cosine);
        assert_float_eq!(*overlap_similarity_already_normalized(&a, &b), *overlap);
    }
}
