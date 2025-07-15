// Modules
mod clustering;
mod misc;
mod simseq;

// Imports
use bio::{io::fastq, pattern_matching::myers::BitVec};
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

    let mut rng = SmallRng::seed_from_u64(2);
    //let mut rng = SmallRng::from_os_rng();

    let myrseqs: Vec<MyrSeq> = [
        simseq::generate_pseudo_amplicon(&mut rng, 759, 210, "1"),
        simseq::generate_pseudo_amplicon(&mut rng, 410, 502, "2"),
    ]
    .concat();

    clustering::method_1(myrseqs)?;

    Ok(())
}
