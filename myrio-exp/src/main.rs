#![allow(unused)]

// Imports
use std::collections::HashMap;

use anyhow::Ok;
use bio_seq::prelude::*;
use myrio_core::{
    data::MyrSeq,
    simseq::{Generator, distr::DiscreteDistribution},
};
use myrio_exp::{
    clustering::partition::compute_cluster_cost, scripts::load_testset, simfunc::SimFunc, tax::basic_test,
};
use rand::{SeedableRng, seq::IndexedRandom};

fn main() -> anyhow::Result<()> {
    myrio_core::tax::store::TaxTreeStore::load_from_fasta_file(
        "./ignore/Magnoliopsida_trnH-psbA_mdb.fasta",
        "trnH-psbA",
        Some(&[5, 10]),
        Some(2),
    )?;
    Ok(())
}

fn cluster_simple_test() {
    let mut rng = rand::rngs::StdRng::from_os_rng();
    let generator = Generator::default()
        .with_q_score_distr(DiscreteDistribution::new_nbin_from_mean_and_std(10.0, 6.0).unwrap());
    let myrseqs: Vec<MyrSeq> = [
        // ITS
        generator.generate_pseudo_amplicon(588, 150, "1", &mut rng),
        // matK
        generator.generate_pseudo_amplicon(1500, 150, "2", &mut rng),
        // rbcL
        generator.generate_pseudo_amplicon(1431, 150, "3", &mut rng),
        // trnH-psbA
        generator.generate_pseudo_amplicon(1058, 150, "4", &mut rng),
    ]
    .concat();
    let clusters =
        myrio_exp::clustering::partition::Clusterer::_cluster_sparse(myrseqs, 4, 4, 0.2, SimFunc::Cosine);
    /*
    for (idx, cluster) in clusters.iter().enumerate() {
        let mut id_count_map: HashMap<&str, usize> = HashMap::new();
        for myrseq in cluster.iter() {
            let count_ref = id_count_map.entry(myrseq.id.as_str()).or_default();
            *count_ref += 1;
        }
        println!("cluster {idx}: {id_count_map:?}")
    }
    */
}

/*
fn simseq_test() -> anyhow::Result<()> {
    let mut rng = rand::rngs::StdRng::from_os_rng();
    let generator = Generator::default();
    let myrseqs = generator.generate_pseudo_amplicon(200, 100, "test", &mut rng);

    for myrseq in myrseqs {
        println!("{}", myrseq.description.unwrap())
    }

    Ok(())
}

fn argmin_run() -> anyhow::Result<()> {
    const K: usize = 5;
    clustering::optim::run_all(K)
}

fn argmin_check() -> anyhow::Result<()> {
    let myrseqs = MyrSeq::decode_vec_from_file("./ignore/argmin_myrseqs.bin")?;

    cluster_and_display(
        myrseqs.clone(),
        unimplemented!(),
        clustering::cosine_similarity,
        5,
        0.6733343002437883,
        0.6041194588723975,
    );

    cluster_and_display(
        myrseqs.clone(),
        unimplemented!(),
        clustering::overlap_similarity,
        5,
        0.6979517680844866,
        0.7561283094321495,
    );

    cluster_and_display(
        myrseqs.clone(),
        clustering::method_two,
        clustering::cosine_similarity,
        5,
        0.6988019466400146,
        0.7003933906555175,
    );

    cluster_and_display(
        myrseqs.clone(),
        clustering::method_two,
        clustering::overlap_similarity,
        5,
        0.3967739229380448,
        0.6880462781710136,
    );

    Ok(())
}

fn cluster_and_display(
    myrseqs: Vec<MyrSeq>,
    cluster_method: ClusterMethod,
    sim_func: SimFunc,
    k: usize,
    t1_cutoff: f64,
    t2_cutoff: f64,
) {
    println!(
        "## cluster_method: {}, sim_func: {}",
        unimplemented!(),
        if std::ptr::fn_addr_eq(
            sim_func,
            clustering::cosine_similarity
                as for<'a, 'b> fn(&'a HashMap<usize, f64>, &'b HashMap<usize, f64>) -> f64
        ) {
            "cosine"
        } else {
            "overlap"
        }
    );
    let clusters = cluster_method(myrseqs.clone(), 5, t1_cutoff, t2_cutoff, sim_func).unwrap();

    for (idx, cluster) in clusters.iter().enumerate() {
        let mut id_count_map: HashMap<&str, usize> = HashMap::new();
        for myrseq in cluster.iter() {
            let count_ref = id_count_map.entry(myrseq.id.as_str()).or_default();
            *count_ref += 1;
        }
        println!("cluster {idx}: {id_count_map:?}")
    }
    println!();
}
*/
