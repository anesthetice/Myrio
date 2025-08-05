use std::{collections::HashMap, ops::AddAssign};

use bio_seq::prelude::*;
use divan::{Bencher, black_box};
use itertools::Itertools;
use myrio_core::data::sparse::SparseFloatVec;
use ndarray::Array1;

const K: usize = 5;
const NB_KMERS: usize = 4_usize.pow(K as u32);

#[divan::bench]
fn kmer_map_array(bencher: Bencher) {
    let seq = black_box(dna!(
        "GCGAGACGAGGGAGTAGGCGCAAGTCGATCGCGCTGTATCCATATATTCATCGTTCCCAAGTAATGTCTCGAAGACATTTTACATAATCGGTCATGCGATGGGAATCGATAGCGGTCAGTGAGCTTAGAGGTCGACTCCAAACGTTAACT"
    ));

    bencher.bench(move || -> [f64; NB_KMERS] {
        let mut map = [0_f64; NB_KMERS];
        for (kmer, kmer_rc) in seq.kmers::<K>().zip_eq(seq.to_revcomp().kmers::<K>()) {
            map[usize::from(&kmer)] += 1.0;
            map[usize::from(&kmer_rc)] += 1.0;
        }
        map
    });
}

#[divan::bench]
fn kmer_map_ndarray(bencher: Bencher) {
    let seq = black_box(dna!(
        "GCGAGACGAGGGAGTAGGCGCAAGTCGATCGCGCTGTATCCATATATTCATCGTTCCCAAGTAATGTCTCGAAGACATTTTACATAATCGGTCATGCGATGGGAATCGATAGCGGTCAGTGAGCTTAGAGGTCGACTCCAAACGTTAACT"
    ));

    bencher.bench(move || -> Array1<f64> {
        let mut map = Array1::<f64>::zeros(NB_KMERS);
        for (kmer, kmer_rc) in seq.kmers::<K>().zip_eq(seq.to_revcomp().kmers::<K>()) {
            map[usize::from(&kmer)] += 1.0;
            map[usize::from(&kmer_rc)] += 1.0;
        }
        map
    });
}

#[divan::bench]
fn kmer_map_vec(bencher: Bencher) {
    let seq = black_box(dna!(
        "GCGAGACGAGGGAGTAGGCGCAAGTCGATCGCGCTGTATCCATATATTCATCGTTCCCAAGTAATGTCTCGAAGACATTTTACATAATCGGTCATGCGATGGGAATCGATAGCGGTCAGTGAGCTTAGAGGTCGACTCCAAACGTTAACT"
    ));

    bencher.bench(move || -> Vec<f64> {
        let mut map: Vec<f64> = vec![0.0; NB_KMERS];
        for (kmer, kmer_rc) in seq.kmers::<K>().zip_eq(seq.to_revcomp().kmers::<K>()) {
            map[usize::from(&kmer)] += 1.0;
            map[usize::from(&kmer_rc)] += 1.0;
        }
        map
    });
}

#[divan::bench]
fn kmer_map_sfvec(bencher: Bencher) {
    let seq = black_box(dna!(
        "GCGAGACGAGGGAGTAGGCGCAAGTCGATCGCGCTGTATCCATATATTCATCGTTCCCAAGTAATGTCTCGAAGACATTTTACATAATCGGTCATGCGATGGGAATCGATAGCGGTCAGTGAGCTTAGAGGTCGACTCCAAACGTTAACT"
    ));

    bencher.bench(move || -> SparseFloatVec {
        let mut map = SparseFloatVec::new();
        for (kmer, kmer_rc) in seq.kmers::<K>().zip_eq(seq.to_revcomp().kmers::<K>()) {
            map.get_or_insert_mut(usize::from(&kmer)).add_assign(1.0);
            map.get_or_insert_mut(usize::from(&kmer_rc)).add_assign(1.0);
        }
        map
    });
}

#[divan::bench]
fn kmer_map_hashmap(bencher: Bencher) {
    let seq = black_box(dna!(
        "GCGAGACGAGGGAGTAGGCGCAAGTCGATCGCGCTGTATCCATATATTCATCGTTCCCAAGTAATGTCTCGAAGACATTTTACATAATCGGTCATGCGATGGGAATCGATAGCGGTCAGTGAGCTTAGAGGTCGACTCCAAACGTTAACT"
    ));

    bencher.bench(move || -> HashMap<usize, f64> {
        let mut map: HashMap<usize, f64> = HashMap::new();
        for (kmer, kmer_rc) in seq.kmers::<K>().zip_eq(seq.to_revcomp().kmers::<K>()) {
            map.entry(usize::from(&kmer)).or_default().add_assign(1.0);
            map.entry(usize::from(&kmer_rc)).or_default().add_assign(1.0);
        }
        map
    });
}

fn main() {
    // Run registered benchmarks.
    divan::main();
}
