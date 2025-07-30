use divan::{Bencher, black_box};
use myrio_core::MyrSeq;

#[divan::bench]
fn hashmap(bencher: Bencher) {}

fn main() {
    // Run registered benchmarks.
    divan::main();
}
