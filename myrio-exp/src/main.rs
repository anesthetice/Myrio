#![allow(unused)]

// Modules
mod clustering;
mod misc;

use std::collections::HashMap;

// Imports
use myrio_core::{MyrSeq, simseq::Generator};
use rand::{SeedableRng, seq::IndexedRandom};

use crate::clustering::{ClusterMethod, SimFunc};

fn main() -> anyhow::Result<()> {
    argmin_check()?;
    Ok(())
}

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
        clustering::method_one,
        clustering::cosine_similarity,
        5,
        0.6733343002437883,
        0.6041194588723975,
    );

    cluster_and_display(
        myrseqs.clone(),
        clustering::method_one,
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
        if std::ptr::fn_addr_eq(
            cluster_method,
            clustering::method_one
                as fn(
                    Vec<MyrSeq>,
                    usize,
                    f64,
                    f64,
                    for<'a, 'b> fn(&'a HashMap<usize, f64>, &'b HashMap<usize, f64>) -> f64,
                ) -> Result<Vec<Vec<MyrSeq>>, anyhow::Error>
        ) {
            "one"
        } else {
            "two"
        },
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
