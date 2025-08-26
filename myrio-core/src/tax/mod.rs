// Modules
pub mod clade;
pub mod store;
pub mod tree;

// Imports
use bio_seq::{
    ReverseComplement,
    codec::{dna::Dna, iupac::Iupac},
    kmer::Kmer,
    seq::{Seq, SeqSlice},
};
use itertools::Itertools;
use myrio_proc::gen_match_k_sparse;

use crate::data::{MyrSeq, SFVec};

pub(crate) const MAX_CONSECUTIVE_N_BEFORE_CUTOFF_DEFAULT: usize = 2;

#[allow(non_snake_case)]
pub fn compute_sparse_kmer_counts_for_fasta_seq(
    seq: &SeqSlice<Iupac>,
    k: usize,
    max_consecutive_N_before_cutoff: usize,
) -> SFVec {
    let mut count_N: usize = 0;
    let expanded = seq
        .into_iter()
        .filter(|nu| {
            if matches!(nu, Iupac::N) {
                count_N += 1;
            } else {
                count_N = 0;
            }
            count_N <= max_consecutive_N_before_cutoff
        })
        .map(|iupac| iupac.to_dna_ext())
        .collect_vec();

    let mut pairs: Vec<(usize, f64)> = Vec::new();

    macro_rules! body {
        ($expanded:expr, $K:expr) => {{
            // A gap (`-`) interrupts the counting
            for group in $expanded.split(|avec| avec.is_empty()) {
                for window in group.windows($K) {
                    let combinations =
                        window.iter().map(|avec| avec.iter()).multi_cartesian_product().collect_vec();
                    let weight: f64 = 1.0 / combinations.len() as f64;

                    combinations.into_iter().for_each(|comb| {
                        let kmer_seq: Seq<Dna> = Seq::from(comb);
                        let kmer: Kmer<Dna, $K> = Kmer::unsafe_from_seqslice(&kmer_seq);
                        pairs.push((usize::from(&kmer), weight));
                        pairs.push((usize::from(&kmer.to_revcomp()), weight))
                    });
                }
            }
        }};
    }
    gen_match_k_sparse!(expanded);
    SFVec::from_unsorted_pairs(pairs, 4_usize.pow(k as u32), 0.0)
}

#[cfg(test)]
mod test {}
