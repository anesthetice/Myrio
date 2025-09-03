#![allow(non_snake_case)]

// Modules
pub mod clade;
pub mod compute;
pub mod core;
pub mod results;
pub mod store;

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
use rand::{
    SeedableRng,
    rngs::{SmallRng, StdRng},
    seq::IndexedRandom,
};
use rayon::{
    iter::{IntoParallelIterator, ParallelBridge, ParallelIterator},
    slice::ParallelSlice,
};
use thiserror::Error;

use crate::{
    data::{Float, MyrSeq, SFVec, SparseVec},
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

pub fn kmer_store_counts_to_kmer_counts(input: (SparseVec<u16>, Float)) -> SFVec {
    let (svec, rescale_factor) = input;
    svec.apply_into(|val| val as Float * rescale_factor)
}

pub fn compute_kmer_counts_for_fasta_seq(
    seq: &SeqSlice<Iupac>,
    k: usize,
    nb_bootstrap_resamples: usize,
    high_pass_filter_flag: bool,
    rng: &mut impl rand::Rng,
) -> SFVec {
    let input =
        compute_kmer_store_counts_for_fasta_seq(seq, k, nb_bootstrap_resamples, high_pass_filter_flag, rng);
    kmer_store_counts_to_kmer_counts(input)
}

pub fn compute_kmer_store_counts_for_fasta_seq(
    seq: &SeqSlice<Iupac>,
    k: usize,
    nb_bootstrap_resamples: usize,
    high_pass_filter_flag: bool,
    rng: &mut impl rand::Rng,
) -> (SparseVec<u16>, Float) {
    fn sample_iupac_nc(
        nc: &Iupac,
        rng: &mut impl rand::Rng,
    ) -> Dna {
        match nc {
            Iupac::A => Dna::A,
            Iupac::C => Dna::C,
            Iupac::G => Dna::G,
            Iupac::T => Dna::T,
            Iupac::R => unsafe { *[Dna::A, Dna::G].choose(rng).unwrap_unchecked() },
            Iupac::Y => unsafe { *[Dna::C, Dna::T].choose(rng).unwrap_unchecked() },
            Iupac::S => unsafe { *[Dna::C, Dna::G].choose(rng).unwrap_unchecked() },
            Iupac::W => unsafe { *[Dna::A, Dna::T].choose(rng).unwrap_unchecked() },
            Iupac::K => unsafe { *[Dna::G, Dna::T].choose(rng).unwrap_unchecked() },
            Iupac::M => unsafe { *[Dna::A, Dna::C].choose(rng).unwrap_unchecked() },
            Iupac::B => unsafe { *[Dna::C, Dna::G, Dna::T].choose(rng).unwrap_unchecked() },
            Iupac::D => unsafe { *[Dna::A, Dna::G, Dna::T].choose(rng).unwrap_unchecked() },
            Iupac::H => unsafe { *[Dna::A, Dna::C, Dna::T].choose(rng).unwrap_unchecked() },
            Iupac::V => unsafe { *[Dna::A, Dna::C, Dna::G].choose(rng).unwrap_unchecked() },
            Iupac::N => unsafe { *[Dna::A, Dna::C, Dna::G, Dna::T].choose(rng).unwrap_unchecked() },
            Iupac::X => unreachable!(),
        }
    }
    let base_nucleotides = seq.iter().filter(|nc| !matches!(nc, Iupac::X)).collect_vec();
    if base_nucleotides.len() < k {
        // Not completely ideal, but trimming empty sfvec entries from the compute tree is even more annoying
        return unsafe { (SparseVec::new_unchecked(vec![0], vec![1], 1, 0), 1.0) };
    }

    let nb_kmers_per_seq = 2 * (base_nucleotides.len() - k + 1);
    let mut singles: Vec<usize> = Vec::with_capacity(nb_bootstrap_resamples * nb_kmers_per_seq);

    macro_rules! body {
        ($seq:expr, $K:expr) => {{
            unsafe {
                let ptr: *mut usize = singles.as_mut_ptr();
                let mut idx: usize = 0;
                for _ in 0..nb_bootstrap_resamples {
                    let seq: Seq<Dna> =
                        Seq::from(&base_nucleotides.iter().map(|nc| sample_iupac_nc(nc, rng)).collect_vec());
                    for key in
                        seq.kmers::<$K>().chain(seq.to_revcomp().kmers::<$K>()).map(|kmer| usize::from(&kmer))
                    {
                        ptr.add(idx).write(key);
                        idx += 1;
                    }
                }
                singles.set_len(idx);
            }
        }};
    }
    gen_match_k_sparse!(_);

    let mut sparse_vec = SparseVec::from_unsorted_singles(singles, 1_u16, 4_usize.pow(k as u32), 0_u16);

    if high_pass_filter_flag {
        let dim = sparse_vec.dim();
        let cutoff = ((nb_bootstrap_resamples + 5) >> 3) as u16; // Equivalent to `(nb_boostraps + 5) / 8` then floor but much more efficient
        let (keys, values): (Vec<usize>, Vec<u16>) =
            sparse_vec.iter_owned().filter(|(_, v)| v > &cutoff).multiunzip();

        sparse_vec = unsafe { SparseVec::new_unchecked(keys, values, dim, 0_u16) }
    }
    let rescale_factor: Float =
        nb_kmers_per_seq as Float / sparse_vec.values().copied().map(u32::from).sum::<u32>() as Float;

    (sparse_vec, rescale_factor)
}

#[allow(non_snake_case)]
pub fn compute_sparse_kmer_counts_for_fasta_seq_old(
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

    let mut pairs: Vec<(usize, Float)> = Vec::new();

    macro_rules! body {
        ($expanded:expr, $K:expr) => {{
            // A gap (`-`) interrupts the counting
            for group in $expanded.split(|avec| avec.is_empty()) {
                for window in group.windows($K) {
                    let combinations =
                        window.iter().map(|avec| avec.iter()).multi_cartesian_product().collect_vec();
                    let weight: Float = 1.0 / combinations.len() as Float;

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
    if pairs.is_empty() {
        // Not completely ideal, but trimming empty sfvec entries from the compute tree is even more annoying
        return unsafe { SFVec::new_unchecked(vec![0], vec![1.0], 1, 0.0) };
    }
    SFVec::from_unsorted_pairs(pairs, 4_usize.pow(k as u32), 0.0)
}

#[cfg(test)]
mod test {}
