// Modules
mod clustering;

use std::hash::{DefaultHasher, Hash, Hasher};

// Imports
use bio::io::fastq;
use bio_seq::prelude::*;
use itertools::Itertools;
use myrio_core::MyrSeq;

const K_VEC_HASHMAP_LIMIT: usize = 11;

fn main() -> anyhow::Result<()> {
    let myrseqs: Vec<MyrSeq> = fastq::Reader::from_file("./ignore/Allium_Ursinum_ITS_barcode82.fastq")?
        .records()
        .filter_map(|record| record.ok().as_ref().map(Into::<MyrSeq>::into))
        .collect();
    println!("{}", myrseqs.len());

    clustering::method_1(myrseqs)
}

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
