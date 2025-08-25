#![allow(unused)]

// Imports
use bio_seq::prelude::*;
use itertools::Itertools;
use myrio_core::data::MyrSeq;
use std::hash::{DefaultHasher, Hash, Hasher};

const K: usize = 5;
const W: usize = 9;

fn get_minimizers(myrseqs: &[MyrSeq]) -> Vec<Vec<(Kmer<Dna, K>, u64)>> {
    myrseqs
        .iter()
        .map(|myrseq| {
            myrseq
                .sequence
                .kmers::<K>()
                .chunks(W)
                .into_iter()
                .map(|chunk| {
                    chunk.into_iter().map(|kmer| (kmer, hash(kmer))).min_by_key(|&(_, hash)| hash).unwrap()
                })
                .collect_vec()
        })
        .collect_vec()
}

fn hash<T: Hash>(seq: T) -> u64 {
    let mut hasher = DefaultHasher::new();
    seq.hash(&mut hasher);
    hasher.finish()
}
