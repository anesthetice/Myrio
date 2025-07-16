// Modules
pub mod testing;

// Imports
use std::{collections::HashMap, ops::AddAssign};

use bio_seq::ReverseComplement;
use itertools::Itertools;
use myrio_core::{MyrSeq, constants::Q_TO_BP_CALL_CORRECT_PROB_MAP};

/// Inspired by isONclust3's method, with quite a few differences, notably we use k-mers instead of minimizers as our expected amplicon length isn't very high (~300-1000 bp). Also we do not neglect `low quality seeds`, we instead assign a quality score to each k-mer that defines its 'strength'.
pub fn method_one(
    myrseqs: Vec<MyrSeq>,
    similarity_function: impl Fn(&HashMap<usize, f64>, &HashMap<usize, f64>) -> f64,
) -> Vec<Vec<MyrSeq>> {
    const K: usize = 5;

    // Step 1
    const T1_CUTOFF: f64 = 0.8;
    let mut myrseqs_extra = myrseqs
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
        .sorted_by(|(.., a), (.., b)| a.total_cmp(b)) // ascending order
        .collect_vec();

    // Step 2
    const T2_CUTOFF: f64 = 0.55;
    struct Cluster {
        reference_seeds: HashMap<usize, f64>,
        elements: Vec<MyrSeq>,
    }
    let mut clusters: Vec<Cluster> = Vec::new();
    {
        let (myrseq, ref_seeds, _) = myrseqs_extra.pop().unwrap();
        clusters.push(Cluster { reference_seeds: ref_seeds, elements: vec![myrseq] });
    }

    while let Some((myrseq, seeds, _)) = myrseqs_extra.pop() {
        match clusters
            .iter_mut()
            .find(|cluster| similarity_function(&seeds, &cluster.reference_seeds) > T2_CUTOFF)
        {
            Some(cluster) => cluster.elements.push(myrseq),
            None => clusters.push(Cluster { reference_seeds: seeds, elements: vec![myrseq] }),
        };
    }

    clusters.into_iter().map(|cluster| cluster.elements).collect_vec()
}

/// More classical approach to clustering, every sequence becomes a cluster and we progressively merge these
pub fn method_two(
    myrseqs: Vec<MyrSeq>,
    similarity_function: impl Fn(&HashMap<usize, f64>, &HashMap<usize, f64>) -> f64,
) -> Vec<Vec<MyrSeq>> {
    const K: usize = 5;

    struct Cluster {
        seeds: HashMap<usize, f64>,
        elements: Vec<MyrSeq>,
    }

    impl Cluster {
        fn merge(
            &mut self,
            other: Self,
        ) {
            self.elements.extend(other.elements);
            for (key, val) in other.seeds.into_iter() {
                let val_ref = self.seeds.entry(key).or_default();
                *val_ref += val;
            }
            self.seeds.iter_mut().for_each(|(_, val)| *val /= 2.0);
        }
    }

    // Step 1
    const T1_CUTOFF: f64 = 0.7;
    let mut clusters: Vec<Cluster> = myrseqs
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
            let mut map: HashMap<usize, f64> = HashMap::new();

            for (idx, (kmer, kmer_rc)) in seq.kmers::<K>().zip_eq(seq_rc.kmers::<K>()).enumerate() {
                // Note: `kmer_rc` is not the reverse complement of `kmer`, it's the `idx`-th k-mer of the reverse complement of the sequence; cleaner implementation if `KmerIter` supported `.rev()` method but it doesn't unfortunately.
                let conf_score = conf_score_per_kmer[idx];
                if conf_score > T1_CUTOFF {
                    let val_ref = map.entry(usize::from(&kmer)).or_default();
                    val_ref.add_assign(conf_score)
                }
                let conf_score = conf_score_per_kmer[nb_kmers - 1 - idx];
                if conf_score > T1_CUTOFF {
                    let val_ref = map.entry(usize::from(&kmer_rc)).or_default();
                    val_ref.add_assign(conf_score)
                }
            }
            Cluster { seeds: map, elements: vec![myrseq] }
        })
        .collect_vec();

    // Step 2
    const T2_CUTOFF: f64 = 0.65;
    let mut halt = false;
    while !halt {
        halt = true;
        let mut i = 0;
        while i < clusters.len() {
            let mut j = i + 1;
            while j < clusters.len() {
                if similarity_function(&clusters[i].seeds, &clusters[j].seeds) > T2_CUTOFF {
                    halt = false;
                    let cluster_to_be_merged = clusters.remove(j);
                    clusters[i].merge(cluster_to_be_merged);
                }
                j += 1;
            }
            i += 1;
        }
    }

    clusters.into_iter().map(|cluster| cluster.elements).collect_vec()
}

pub fn cosine_similarity(
    a: &HashMap<usize, f64>,
    b: &HashMap<usize, f64>,
) -> f64 {
    let mut dot = 0.0;
    let mut norm_a = 0.0;
    let mut norm_b = 0.0;

    for (&key, &val_a) in a {
        norm_a += val_a * val_a;
        if let Some(&val_b) = b.get(&key) {
            dot += val_a * val_b;
        }
    }

    for &val_b in b.values() {
        norm_b += val_b * val_b;
    }

    if norm_a == 0.0 || norm_b == 0.0 {
        return 0.0;
    }

    dot / (norm_a.sqrt() * norm_b.sqrt())
}

pub fn overlap_similarity(
    seq: &HashMap<usize, f64>,
    ref_seq: &HashMap<usize, f64>,
) -> f64 {
    let mut seq_total_sum = 0.0;
    let mut seq_overlap_sum = 0.0;

    for (&key, &seq_val) in seq {
        seq_total_sum += seq_val;
        if let Some(&ref_seq_val) = ref_seq.get(&key) {
            seq_overlap_sum += f64::min(seq_val, ref_seq_val)
        }
    }
    seq_overlap_sum / seq_total_sum
}
