// Modules
mod clustering;
mod misc;
mod simseq;

use std::collections::HashMap;

// Imports
use bio::{io::fastq, pattern_matching::myers::BitVec};
use bio_seq::prelude::*;
use itertools::Itertools;
use myrio_core::MyrSeq;
use rand::{RngCore, SeedableRng, rngs::SmallRng};
use rand_distr::{Distribution, Poisson};

fn main() -> anyhow::Result<()> {
    /*
    let myrseqs: Vec<MyrSeq> = myrio_core::io::read_fastq(
        //"./ignore/Cymbopogon_Citrus_Qiagen_matk_rbcL_psbA-trnH_ITS_barcode3.fastq",
        "./ignore/Begonia_Masoniana_Qiagen_matk_rbcL_psbA-trnH_ITS_barcode8.fastq",
    )?;
    println!("→ collected: {} sequences\n", myrseqs.len());
    */

    let mut rng = SmallRng::seed_from_u64(0);

    let myrseqs: Vec<MyrSeq> = [
        simseq::generate_pseudo_amplicon(&mut rng, 759, 210, "1"),
        simseq::generate_pseudo_amplicon(&mut rng, 239, 502, "2"),
        simseq::generate_pseudo_amplicon(&mut rng, 759, 50, "3"),
        simseq::generate_pseudo_amplicon(&mut rng, 902, 50, "4"),
    ]
    .concat();

    let clusters = clustering::method_one(myrseqs.clone(), 5, clustering::cosine_similarity)?;

    for (idx, cluster) in clusters.iter().enumerate() {
        ///*
        let mut id_count_map: HashMap<&str, usize> = HashMap::new();
        for myrseq in cluster.iter() {
            let count_ref = id_count_map.entry(myrseq.id.as_str()).or_default();
            *count_ref += 1;
        }
        println!("cluster {idx}: {id_count_map:?}")
        //*/
        //println!("cluster {idx}: {}", cluster.len());
    }

    std::thread::sleep(std::time::Duration::from_micros(3000));
    println!("\n\n");

    let clusters = clustering::method_two(myrseqs, 5, clustering::cosine_similarity)?;

    for (idx, cluster) in clusters.iter().enumerate() {
        //*/
        let mut id_count_map: HashMap<&str, usize> = HashMap::new();
        for myrseq in cluster.iter() {
            let count_ref = id_count_map.entry(myrseq.id.as_str()).or_default();
            *count_ref += 1;
        }
        println!("cluster {idx}: {id_count_map:?}")
        //*/
        //println!("cluster {idx}: {}", cluster.len());
    }

    Ok(())
}
