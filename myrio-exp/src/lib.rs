// Modules
pub mod clustering;
pub mod misc;
pub mod scripts;
pub mod simfunc;
pub mod simseq;
pub mod tax;

// Imports
use bio_seq::prelude::*;
use itertools::Itertools;
use myrio_core::data::{Float, MyrSeq};
use myrio_proc::gen_match_k_dense;

pub type DFArray = ndarray::Array1<f64>;

pub const K_DENSE_VALID_RANGE: std::ops::Range<usize> = 2..7;
pub const K_DENSE_VALID_RANGE_ERROR_MSG: &str =
    "For dense k-mer count maps, only k ∈ {{2, ..., 6}} is currently supported";

pub fn compute_dense_kmer_counts(
    myrseq: &MyrSeq,
    k: usize,
    cutoff: Float,
) -> (DFArray, usize) {
    let conf_score_per_kmer = myrseq
        .quality
        .iter()
        .map(|q| myrio_core::constants::Q_TO_BP_CALL_CORRECT_PROB_MAP[*q as usize])
        .collect_vec()
        .windows(k)
        .map(|vals| vals.iter().product::<Float>())
        .collect_vec();

    macro_rules! body {
        ($seq:expr, $K:expr) => {{
            let nb_kmers = $seq.len() - $K + 1;
            let mut nb_hck: usize = 0; // number of high-quality k-mers
            let mut map = DFArray::zeros(4_usize.pow(k as u32));

            for (idx, (kmer, kmer_rc)) in $seq.kmers::<$K>().zip_eq($seq.to_revcomp().kmers::<$K>()).enumerate() {
                // Note: `kmer_rc` is not the reverse complement of `kmer`, it's the `idx`-th k-mer of the reverse complement of the sequence.
                if conf_score_per_kmer[idx] > cutoff {
                    map[usize::from(&kmer)] += 1.0;
                    nb_hck += 1;
                }
                if conf_score_per_kmer[nb_kmers - 1 - idx] > cutoff {
                    map[usize::from(&kmer_rc)] += 1.0;
                    nb_hck += 1;
                }
            }
            (map, nb_hck)
        }};
    }
    gen_match_k_dense!(myrseq.sequence)
}
