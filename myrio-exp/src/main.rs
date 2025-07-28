// Modules
mod clustering;
mod misc;

use std::collections::HashMap;

// Imports
use myrio_core::MyrSeq;

fn main3() {
    clustering::testing::run().unwrap();
}

fn main() -> anyhow::Result<()> {
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
