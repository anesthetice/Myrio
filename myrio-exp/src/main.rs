// Modules
mod clustering;
mod misc;

use std::collections::HashMap;

// Imports
use myrio_core::{MyrSeq, simseq::Generator};
use rand::{SeedableRng, seq::IndexedRandom};

fn main() -> anyhow::Result<()> {
    simseq_test()?;
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
    clustering::testing::run()
}

fn argmin_check() -> anyhow::Result<()> {
    let myrseqs = MyrSeq::decode_vec_from_file("./ignore/argmin_myrseqs.bin")?;
    let t1_cutoff = 0.6412506471605364;
    let t2_cutoff = 0.6449838846650066;

    let clusters =
        clustering::method_one(myrseqs.clone(), 5, t1_cutoff, t2_cutoff, clustering::cosine_similarity)?;

    for (idx, cluster) in clusters.iter().enumerate() {
        // /*
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
