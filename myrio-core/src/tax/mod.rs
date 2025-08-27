#![allow(non_snake_case)]

// Modules
pub mod clade;
pub mod core;
pub mod store;
pub mod tree;

// Imports
use bio_seq::{
    ReverseComplement,
    codec::{dna::Dna, iupac::Iupac},
    error::ParseBioError,
    kmer::Kmer,
    seq::{Seq, SeqSlice},
};
use itertools::Itertools;
use myrio_proc::gen_match_k_sparse;
use thiserror::Error;

use crate::{
    data::{MyrSeq, SFVec},
    tax::clade::Rank,
};

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    BincodeDecode(#[from] bincode::error::DecodeError),
    #[error(transparent)]
    BincodeEncode(#[from] bincode::error::EncodeError),
    #[error(transparent)]
    IO(#[from] std::io::Error),
    #[error("Bio-seq parsing error for the record starting on line {1}, {0}")]
    BioSeq(ParseBioError, usize),
    #[error("Failed to parse taxonomic identity of the record starting on line {1}, {0}")]
    CladeParse(clade::ParsingError, usize),
    #[error("Expected the highest rank to be {0}, got {1} instead for the record starting on line {2}")]
    HighestRankMismatch(Rank, Rank, usize),
}

#[allow(non_snake_case)]
pub fn compute_sparse_kmer_counts_for_fasta_seq(
    seq: &SeqSlice<Iupac>,
    k: usize,
    max_consecutive_N_before_gap: usize,
) -> SFVec {
    let mut count_N: usize = 0;
    let mut gap_already_inserted_flag: bool = false;
    let expanded = seq
        .into_iter()
        .filter_map(|nu| {
            if matches!(nu, Iupac::N) {
                count_N += 1;
            } else {
                count_N = 0;
            }
            if count_N >= max_consecutive_N_before_gap {
                if !gap_already_inserted_flag {
                    gap_already_inserted_flag = true;
                    Some(Iupac::X)
                } else {
                    None
                }
            } else {
                gap_already_inserted_flag = false;
                Some(nu)
            }
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
