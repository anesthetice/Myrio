// Modules
mod clustering;
mod misc;
mod simseq;

// Imports
use bio::io::fastq;
use bio_seq::prelude::*;
use itertools::Itertools;
use myrio_core::MyrSeq;
use rand::{RngCore, SeedableRng, rngs::SmallRng};
use rand_distr::{Distribution, Poisson};

fn main() -> anyhow::Result<()> {
    /*
    let myrseqs: Vec<MyrSeq> =
        fastq::Reader::from_file("./ignore/Cymbopogon_Citrus_Qiagen_matk_rbcL_psbA-trnH_ITS_barcode3.fastq")?
            .records()
            .filter_map(|record| record.ok().as_ref().map(Into::<MyrSeq>::into))
            .collect();
    println!("{}", myrseqs.len());

    clustering::method_1(myrseqs)
    */
    //let mut rng = SmallRng::seed_from_u64(0);
    //simseq::generate_pseudo_amplicon(&mut rng, 300, 100);
    println!("{}", 2.5_f64.round() as usize);
    Ok(())
}
