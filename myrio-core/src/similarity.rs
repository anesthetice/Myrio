// Imports
use itertools::Itertools;
use nutype::nutype;

use crate::data::{Float, SFVec};

pub type SimFunc = fn(&SFVec, &SFVec) -> SimScore;

#[nutype(
    default = 0_f32,
    new_unchecked,
    validate(finite),
    derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deref, TryFrom, Default)
)]
pub struct SimScore(f32);

impl std::fmt::Display for SimScore {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "{:.3}", **self)
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
    (a.dot(b) / (a.norm_l2() * b.norm_l2())).try_into().unwrap_or_default()
}

#[cfg(debug_assertions)]
fn cosine(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    (a.dot(b) / (a.norm_l2() * b.norm_l2()))
        .try_into()
        .inspect_err(|e| {
            eprintln!(
                "{e} a.count={} ({} NaN), b.count={} ({} NaN)",
                a.count(),
                a.values().filter(|v| v.is_nan()).count(),
                b.count(),
                b.values().filter(|v| v.is_nan()).count(),
            )
        })
        .unwrap_or_default()
}

#[cfg(not(debug_assertions))]
fn cosine_pre_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    unsafe { SimScore::new_unchecked(a.dot(b)) }
}

#[cfg(debug_assertions)]
fn cosine_pre_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    a.dot(b)
        .try_into()
        .inspect_err(|e| {
            eprintln!(
                "{e} a.count={} ({} NaN), b.count={} ({} NaN)",
                a.count(),
                a.values().filter(|v| v.is_nan()).count(),
                b.count(),
                b.values().filter(|v| v.is_nan()).count(),
            )
        })
        .unwrap_or_default()
}

#[cfg(not(debug_assertions))]
fn jacard_tanimoto(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.dot(b);
    (dot / (a.norm_l2_squared() + b.norm_l2_squared() - dot)).try_into().unwrap_or_default()
}

#[cfg(debug_assertions)]
fn jacard_tanimoto(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.dot(b);
    (dot / (a.norm_l2_squared() + b.norm_l2_squared() - dot))
        .try_into()
        .inspect_err(|e| {
            eprintln!(
                "{e} a.count={} ({} NaN), b.count={} ({} NaN)",
                a.count(),
                a.values().filter(|v| v.is_nan()).count(),
                b.count(),
                b.values().filter(|v| v.is_nan()).count(),
            )
        })
        .unwrap_or_default()
}

#[cfg(not(debug_assertions))]
fn jacard_tanimoto_pre_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.dot(b);
    (dot / (2.0 - dot)).try_into().unwrap_or_default()
}

#[cfg(debug_assertions)]
fn jacard_tanimoto_pre_normalized(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let dot = a.dot(b);
    (dot / (2.0 - dot))
        .try_into()
        .inspect_err(|e| {
            eprintln!(
                "{e} a.count={} ({} NaN), b.count={} ({} NaN)",
                a.count(),
                a.values().filter(|v| v.is_nan()).count(),
                b.count(),
                b.values().filter(|v| v.is_nan()).count(),
            )
        })
        .unwrap_or_default()
}

#[cfg(not(debug_assertions))]
fn overlap(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let min_max = a.merge_and_apply(b, |x, y| (x.min(y), x.max(y)));

    let sum_min: Float = min_max.values().map(|(min, _)| min).sum();
    let sum_max: Float = min_max.values().map(|(_, max)| max).sum();

    (sum_min / sum_max).try_into().unwrap_or_default()
}

#[cfg(debug_assertions)]
fn overlap(
    a: &SFVec,
    b: &SFVec,
) -> SimScore {
    let min_max = a.merge_and_apply(b, |x, y| (x.min(y), x.max(y)));

    let sum_min: Float = min_max.values().map(|(min, _)| min).sum();
    let sum_max: Float = min_max.values().map(|(_, max)| max).sum();

    (sum_min / sum_max)
        .try_into()
        .inspect_err(|e| {
            eprintln!(
                "{e} a.count={} ({} NaN), b.count={} ({} NaN)",
                a.count(),
                a.values().filter(|v| v.is_nan()).count(),
                b.count(),
                b.values().filter(|v| v.is_nan()).count(),
            )
        })
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

        let mut a = myrseq_1.compute_kmer_counts(2, Float::MIN).0;
        let mut b = myrseq_2.compute_kmer_counts(2, Float::MIN).0;
        let cosine = cosine(&a, &b);
        assert_float_eq!(*cosine, 2.0 / 3.0);
        let overlap = overlap(&a, &b);
        assert_float_eq!(*overlap, 0.5);

        a.normalize_l2();
        b.normalize_l2();

        assert_float_eq!(*cosine_pre_normalized(&a, &b), *cosine);
    }
}
