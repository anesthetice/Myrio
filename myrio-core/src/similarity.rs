// Imports
use itertools::Itertools;
use nutype::nutype;

use crate::data::{Float, SFVec};

pub type SimFunc = fn(&SFVec, &SFVec) -> SimScore;

#[nutype(
    default = 0_f32,
    new_unchecked,
    validate(finite),
    derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deref, TryFrom, Display, Default)
)]
pub struct SimScore(f32);

#[derive(Debug, Clone, Copy)]
pub enum Similarity {
    Cosine,
    /// Called both the `Jaccard index` and the `Tanimoto coefficient`
    /// See `https://link.springer.com/article/10.1007/s10910-010-9668-4`for more info
    JacardTanimoto,
    Overlap,
}

impl Similarity {
    pub fn normalization_improves_performance(&self) -> bool {
        matches!(self, Self::Cosine | Self::JacardTanimoto)
    }

    pub fn to_simfunc(
        &self,
        is_input_normalized: bool,
    ) -> SimFunc {
        match (self, is_input_normalized) {
            (Self::Cosine, true) => cosine_pre_normalized,
            (Self::Cosine, false) => cosine,
            (Self::JacardTanimoto, true) => jacard_tanimoto_pre_normalized,
            (Self::JacardTanimoto, false) => jacard_tanimoto,
            (Self::Overlap, _) => overlap,
        }
    }

    pub fn dist(score: SimScore) -> SimScore {
        // might be necessary to `match self` if other similarity methods are added, currently isn't
        unsafe { SimScore::new_unchecked(1.0 - *score) }
    }
}

#[cfg(not(debug_assertions))]
fn cosine(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.merge_and_apply(b, |x, y| x * y).sum();
    (dot / (a.norm_l2() * b.norm_l2())).try_into().unwrap_or_default()
}

#[cfg(debug_assertions)]
fn cosine(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.merge_and_apply(b, |x, y| x * y).sum();
    (dot / (a.norm_l2() * b.norm_l2()))
        .try_into()
        .inspect_err(|e| eprintln!("{e}, a.count={}, b.count={}", a.count(), b.count()))
        .unwrap_or_default()
}

#[cfg(not(debug_assertions))]
fn cosine_pre_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    unsafe { SimScore::new_unchecked(a.merge_and_apply(b, |x, y| x * y).sum()) }
}

#[cfg(debug_assertions)]
fn cosine_pre_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    a.merge_and_apply(b, |x, y| x * y)
        .sum()
        .try_into()
        .inspect_err(|e| eprintln!("{e}, a.count={}, b.count={}", a.count(), b.count()))
        .unwrap_or_default()
}

#[cfg(not(debug_assertions))]
fn jacard_tanimoto(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.merge_and_apply(b, |x, y| x * y).sum();
    (dot / (a.norm_l2_squared() + b.norm_l2_squared() - dot)).try_into().unwrap_or_default()
}

#[cfg(debug_assertions)]
fn jacard_tanimoto(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.merge_and_apply(b, |x, y| x * y).sum();
    (dot / (a.norm_l2_squared() + b.norm_l2_squared() - dot))
        .try_into()
        .inspect_err(|e| eprintln!("{e}, a.count={}, b.count={}", a.count(), b.count()))
        .unwrap_or_default()
}

#[cfg(not(debug_assertions))]
fn jacard_tanimoto_pre_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.merge_and_apply(b, |x, y| x * y).sum();
    (dot / (2.0 - dot)).try_into().unwrap_or_default()
}

#[cfg(debug_assertions)]
fn jacard_tanimoto_pre_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.merge_and_apply(b, |x, y| x * y).sum();
    (dot / (2.0 - dot))
        .try_into()
        .inspect_err(|e| eprintln!("{e}, a.count={}, b.count={}", a.count(), b.count()))
        .unwrap_or_default()
}

#[cfg(not(debug_assertions))]
fn overlap(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let min_max = a.merge_and_apply(b, |x, y| (x.min(y), x.max(y)));
    let (min_s, max_s) = min_max.sval();

    let nb_sparse = min_max.count() as Float;
    let sum_min = min_max.values().map(|(min, _)| min).sum1().unwrap_or(0.0) + min_s * nb_sparse;
    let sum_max = min_max.values().map(|(_, max)| max).sum1().unwrap_or(0.0) + max_s * nb_sparse;

    (sum_min / sum_max).try_into().unwrap_or_default()
}

#[cfg(debug_assertions)]
fn overlap(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let min_max = a.merge_and_apply(b, |x, y| (x.min(y), x.max(y)));
    let (min_s, max_s) = min_max.sval();

    let nb_sparse = min_max.count() as Float;
    let sum_min = min_max.values().map(|(min, _)| min).sum1().unwrap_or(0.0) + min_s * nb_sparse;
    let sum_max = min_max.values().map(|(_, max)| max).sum1().unwrap_or(0.0) + max_s * nb_sparse;

    (sum_min / sum_max)
        .try_into()
        .inspect_err(|e| eprintln!("{e}, a.count={}, b.count={}", a.count(), b.count()))
        .unwrap_or_default()
}

#[cfg(test)]
mod test {

    use bio_seq::prelude::*;

    use super::*;
    use crate::{assert_float_eq, data::MyrSeq};

    #[test]
    fn simfunc_test() {
        let myrseq_1 = MyrSeq::create("", None, dna!["ACTG"], &[0; 4]);
        let myrseq_2 = MyrSeq::create("", None, dna!["ACTC"], &[0; 4]);

        let mut a = myrseq_1.compute_sparse_kmer_counts(2, f64::MIN).0;
        let mut b = myrseq_2.compute_sparse_kmer_counts(2, f64::MIN).0;
        let cosine = cosine(&a, &b);
        assert_float_eq!(*cosine, 2.0 / 3.0);
        let overlap = overlap(&a, &b);
        assert_float_eq!(*overlap, 0.5);

        a.normalize_l2();
        b.normalize_l2();

        assert_float_eq!(*cosine_pre_normalized(&a, &b), *cosine);
    }
}
