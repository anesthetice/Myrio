use std::{
    io::{BufRead, BufReader},
    path::Path,
    str::FromStr,
};

use anyhow::Context;
use bio_seq::{codec::dna::Dna as DnaCodec, seq::Seq};
use itertools::Itertools;
use myrio_core::{clustering::SimFunc, data::MyrSeq};

pub fn basic_test() -> anyhow::Result<()> {
    const K: usize = 8;

    let refs = std::fs::read_to_string("./ignore/Magnoliopsida_rbcL_raxdb.fasta")?
        .lines()
        .tuple_windows::<(_, _)>()
        .enumerate()
        .step_by(2)
        .map(|(line_idx, (l1, l2))| {
            if !l1.starts_with(">") {
                panic!();
            }
            let (id, desc) = match l1[1..].split_once("|") {
                Some((id, desc)) => (id.to_string(), Some(desc.to_string())),
                None => (l1[1..].to_string(), None),
            };
            let l2: String =
                l2.chars().filter(|c| *c == 'A' || *c == 'C' || *c == 'T' || *c == 'G').collect();
            let seq: Seq<DnaCodec> = Seq::from_str(&l2).with_context(|| line_idx.to_string()).unwrap();
            let len = seq.len();
            let mseq = MyrSeq::new(id, desc, seq, vec![1; len]);
            let kmer_counts = mseq.compute_dense_kmer_counts(K, 0.0).unwrap().0;
            (mseq.id, mseq.description, kmer_counts)
        })
        .collect_vec();

    let queries =
        myrio_core::io::read_fastq_from_file("./ignore/Solanum_lycopersicummatK_rbcL_ITS_barcode11.fastq")?;

    let centroids =
        crate::clustering::partition::Clusterer::cluster_dense_alt(queries, 2, K, 0.2, SimFunc::Cosine);

    for centroid in centroids.into_iter() {
        let best = refs
            .iter()
            .max_by_key(|(_, _, kcount)| SimFunc::Overlap.compute_dense(&centroid, kcount))
            .unwrap();

        println!(
            "→ {}  with score of {:.4}",
            best.1.as_ref().unwrap(),
            SimFunc::Overlap.compute_dense(&centroid, &best.2)
        )
    }

    Ok(())
}
