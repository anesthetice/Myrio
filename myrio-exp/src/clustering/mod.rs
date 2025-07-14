use std::{collections::HashMap, ops::AddAssign};

use bio_seq::ReverseComplement;
use itertools::Itertools;
use myrio_core::{MyrSeq, constants::Q_TO_BP_CALL_CORRECT_PROB_MAP};

pub fn method_1(myrseqs: Vec<MyrSeq>) -> anyhow::Result<()> {
    const K: usize = 5;

    // step 1:
    const T1_CUTOFF: f64 = 0.8;
    let myrseqs_extra = myrseqs
        .into_iter()
        .map(|myrseq| {
            let seq = &myrseq.sequence;
            let seq_rc = seq.to_revcomp();

            let conf_score_per_kmer = myrseq
                .quality
                .iter()
                .map(|q| Q_TO_BP_CALL_CORRECT_PROB_MAP[*q as usize])
                .collect_vec()
                .windows(K)
                .map(|vals| vals.iter().product::<f64>())
                .collect_vec();

            let nb_kmers = seq.len() - K + 1;

            #[allow(non_snake_case)]
            let mut nb_HCS_weighted: f64 = 0.0;
            let mut map: HashMap<usize, f64> = HashMap::new();

            for (idx, (kmer, kmer_rc)) in seq.kmers::<K>().zip_eq(seq_rc.kmers::<K>()).enumerate() {
                // Note: `kmer_rc` is not the reverse complement of `kmer`, it's the `idx`-th k-mer of the reverse complement of the sequence; cleaner implementation if `KmerIter` supported `.rev()` method but it doesn't unfortunately.
                let conf_score = conf_score_per_kmer[idx];
                if conf_score > T1_CUTOFF {
                    let val_ref = map.entry(usize::from(&kmer)).or_default();
                    nb_HCS_weighted += conf_score;
                    val_ref.add_assign(conf_score)
                }
                let conf_score = conf_score_per_kmer[nb_kmers - 1 - idx];
                if conf_score > T1_CUTOFF {
                    let val_ref = map.entry(usize::from(&kmer_rc)).or_default();
                    nb_HCS_weighted += conf_score;
                    val_ref.add_assign(conf_score)
                }
            }
            (myrseq, map, nb_HCS_weighted)
        })
        .sorted_by(|(.., a), (.., b)| b.total_cmp(a))
        .collect_vec();

    for (.., c) in a {
        println!("{c}")
    }

    Ok(())
}
