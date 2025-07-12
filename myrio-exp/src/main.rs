// Modules
mod clustering;

use std::{
    hash::{DefaultHasher, Hash, Hasher},
    path::Path,
};

// Imports
use bio::io::fastq;
use bio_seq::{codec::Codec, prelude::*};
use itertools::Itertools;
use myrio_core::MyrSeq;

const K_VEC_HASHMAP_LIMIT: usize = 11;

const K: usize = 4;
const W: usize = 9;

fn main() {
    #[allow(arithmetic_overflow)]
    let a = 1_usize << (K * bio_seq::codec::dna::Dna::BITS as usize);
    println!("{}", a);

    let b = dna!("ACGT");
    b.kmers::<2>().for_each(|kmer| println!("{kmer:?}"));
}

fn main2() -> anyhow::Result<()> {
    let myrseqs: Vec<MyrSeq> = fastq::Reader::from_file("./ignore/Allium_Ursinum_ITS_barcode82.fastq")?
        .records()
        .filter_map(|record| record.ok().as_ref().map(Into::<MyrSeq>::into))
        .collect();
    println!("{}", myrseqs.len());

    let super_minimizers = myrseqs
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
        .collect_vec();

    let myrseq0 = myrseqs[0].clone();
    let myrseq0_minimizers = super_minimizers[0].clone();
    println!("{:#?}\n\n{}", myrseq0, myrseq0_minimizers.iter().map(|a| a.0.to_string()).join(", "));

    Ok(())
}

fn hash<T: Hash>(seq: T) -> u64 {
    let mut hasher = DefaultHasher::new();
    seq.hash(&mut hasher);
    hasher.finish()
}
