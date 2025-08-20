// Imports
use divan::{Bencher, black_box};
use myrio_core::{
    clustering::{SimFunc, SimScore},
    data::{DFArray, MyrSeq, SFVec},
    simseq::Generator,
};
use rand::SeedableRng;

const K: usize = 6;
const SEED: u64 = 0;
const LENGTH: usize = 500;

#[divan::bench]
fn dense_cosine(bencher: Bencher) {
    let mut rng = rand::rngs::SmallRng::seed_from_u64(SEED);
    let a: DFArray = black_box(
        MyrSeq::create("a", None, &Generator::generate_core_sequence(LENGTH, &mut rng), &vec![1; LENGTH])
            .compute_dense_kmer_counts(K, 0.0)
            .unwrap()
            .0,
    );
    let b: DFArray = black_box(
        MyrSeq::create("b", None, &Generator::generate_core_sequence(LENGTH, &mut rng), &vec![1; LENGTH])
            .compute_dense_kmer_counts(K, 0.0)
            .unwrap()
            .0,
    );

    bencher.bench(move || -> SimScore { SimFunc::Cosine.compute_dense(&a, &b) });
}

#[divan::bench]
fn sparse_cosine(bencher: Bencher) {
    let mut rng = rand::rngs::SmallRng::seed_from_u64(SEED);
    let a: SFVec = black_box(
        MyrSeq::create("a", None, &Generator::generate_core_sequence(LENGTH, &mut rng), &vec![1; LENGTH])
            .compute_sparse_kmer_counts(K, 0.0)
            .unwrap()
            .0,
    );
    let b: SFVec = black_box(
        MyrSeq::create("b", None, &Generator::generate_core_sequence(LENGTH, &mut rng), &vec![1; LENGTH])
            .compute_sparse_kmer_counts(K, 0.0)
            .unwrap()
            .0,
    );

    bencher.bench(move || SimFunc::Cosine.compute_sparse(&a, &b));
}

#[divan::bench]
fn dense_overlap(bencher: Bencher) {
    let mut rng = rand::rngs::SmallRng::seed_from_u64(SEED);
    let a: DFArray = black_box(
        MyrSeq::create("a", None, &Generator::generate_core_sequence(LENGTH, &mut rng), &vec![1; LENGTH])
            .compute_dense_kmer_counts(K, 0.0)
            .unwrap()
            .0,
    );
    let b: DFArray = black_box(
        MyrSeq::create("b", None, &Generator::generate_core_sequence(LENGTH, &mut rng), &vec![1; LENGTH])
            .compute_dense_kmer_counts(K, 0.0)
            .unwrap()
            .0,
    );

    bencher.bench(move || SimFunc::Overlap.compute_dense(&a, &b));
}

#[divan::bench]
fn sparse_overlap(bencher: Bencher) {
    let mut rng = rand::rngs::SmallRng::seed_from_u64(SEED);
    let a: SFVec = black_box(
        MyrSeq::create("a", None, &Generator::generate_core_sequence(LENGTH, &mut rng), &vec![1; LENGTH])
            .compute_sparse_kmer_counts(K, 0.0)
            .unwrap()
            .0,
    );
    let b: SFVec = black_box(
        MyrSeq::create("b", None, &Generator::generate_core_sequence(LENGTH, &mut rng), &vec![1; LENGTH])
            .compute_sparse_kmer_counts(K, 0.0)
            .unwrap()
            .0,
    );

    bencher.bench(move || SimFunc::Overlap.compute_sparse(&a, &b));
}
