use itertools::Itertools;
use rayon::{iter::ParallelIterator, slice::ParallelSlice};

use crate::{
    data::{MyrSeq, SFVec},
    similarity::{self, SimFunc, SimScore},
};

pub struct Clusterer {
    myrseqs: Vec<MyrSeq>,
    k: usize,
    t1: f64,
    t2: usize,
}

impl Clusterer {
    pub fn new(
        myrseqs: Vec<MyrSeq>,
        k: usize,
        t1: f64,
        t2: usize,
    ) -> Self {
        Self { myrseqs, k, t1, t2 }
    }

    pub fn compute_cluster_centroid<T>(elements: &[T]) -> SFVec
    where
        T: AsRef<SFVec> + Sync,
    {
        let n = elements.len() as f64;

        let par_chunks = elements.par_chunks_exact(2);
        let remainder = par_chunks.remainder().first();

        // We have to use a slightly annoying trick as reduce_with expects same input output type
        let sfvec = par_chunks
            .map(|chunk| unsafe {
                chunk
                    .get_unchecked(0)
                    .as_ref()
                    .merge_with_and_apply(chunk.get_unchecked(1).as_ref(), |x, y| x + y)
            })
            .reduce_with(|a, b| a.merge_with_and_apply(&b, |x, y| x + y))
            .unwrap();

        if let Some(rem) = remainder {
            sfvec.merge_with_and_apply(rem.as_ref(), |x, y| x + y) / n
        } else {
            sfvec / n
        }
    }

    pub fn cluster_with_initial_number_of_cluster(
        self,
        nb_clusters: usize,
    ) -> Vec<Vec<MyrSeq>> {
        let simfunc = similarity::cosine_similarity_already_normalized;

        // Step 1, compute the k-mer freqs of each myrseq
        let (myrseqs, kmer_freqs_vec, nb_hck_vec): (Vec<MyrSeq>, Vec<SFVec>, Vec<usize>) = self
            .myrseqs
            .into_iter()
            .filter_map(|myrseq| {
                let (mut kmer_freqs, nb_hck) = myrseq.compute_sparse_kmer_counts(self.k, self.t1);
                if nb_hck > self.t2 {
                    kmer_freqs /= nb_hck as f64;
                    Some((myrseq, kmer_freqs, nb_hck))
                } else {
                    None
                }
            })
            //.sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
            .multiunzip();

        // Step 2, initialize the clusters
        let mut clusters: Vec<Cluster> = vec![Cluster::new(&kmer_freqs_vec[0])];
        let max_nb_hck = nb_hck_vec[0] as f64;
        while clusters.len() < nb_clusters {
            let new_cluster_seeds = kmer_freqs_vec
                .iter()
                .zip_eq(nb_hck_vec.iter())
                .min_by_key(|&(kcount, nb_hck)| {
                    let mut score: f64 =
                        clusters.iter().map(|cluster| *simfunc(&cluster.centroid_seeds, kcount)).sum();

                    // mean similarity score
                    score /= clusters.len() as f64;
                    // penalty (increases score) for seqs with a low number of high-confidence k-mers
                    score = 0.7 * score + 0.3 * (1.0 - (*nb_hck as f64) / max_nb_hck);
                    SimScore::try_new(score).unwrap()
                })
                .unwrap()
                .0;
            clusters.push(Cluster::new(new_cluster_seeds));
        }

        Self::cluster_core(myrseqs, kmer_freqs_vec, clusters)
    }

    pub fn cluster_with_initial_centroids(
        self,
        initial_centroids: Vec<SFVec>,
    ) -> Vec<Vec<MyrSeq>> {
        // Step 1, compute the k-mer freqs of each myrseq
        let (myrseqs, kmer_freqs_vec, nb_hck_vec): (Vec<MyrSeq>, Vec<SFVec>, Vec<usize>) = self
            .myrseqs
            .into_iter()
            .filter_map(|myrseq| {
                let (mut kmer_freqs, nb_hck) = myrseq.compute_sparse_kmer_counts(self.k, self.t1);
                if nb_hck > self.t2 {
                    kmer_freqs /= nb_hck as f64;
                    Some((myrseq, kmer_freqs, nb_hck))
                } else {
                    None
                }
            })
            //.sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
            .multiunzip();

        // Step 2, initialize the clusters
        let clusters = initial_centroids
            .into_iter()
            .map(|centroid| Cluster { centroid_seeds: centroid, elements: Vec::new() })
            .collect_vec();

        Self::cluster_core(myrseqs, kmer_freqs_vec, clusters)
    }

    fn cluster_core(
        myrseqs: Vec<MyrSeq>,
        kmer_freqs_vec: Vec<SFVec>,
        initial_clusters: Vec<Cluster>,
    ) -> Vec<Vec<MyrSeq>> {
        let simfunc = similarity::cosine_similarity_already_normalized;
        let mut clusters = initial_clusters;
        // Step 3
        for _ in 0..3 {
            for kmer_freqs in kmer_freqs_vec.iter() {
                clusters
                    .iter_mut()
                    .max_by_key(|cluster| simfunc(&cluster.centroid_seeds, kmer_freqs))
                    .unwrap()
                    .elements
                    .push(kmer_freqs);
            }

            clusters = clusters.into_iter().map(|cluster| cluster.recenter()).collect_vec();
        }

        // Step 4
        let mut output: Vec<Vec<MyrSeq>> = vec![Vec::new(); clusters.len()];
        for (myrseq, kmer_freqs) in myrseqs.into_iter().zip_eq(kmer_freqs_vec.iter()) {
            let idx = clusters
                .iter()
                .enumerate()
                .max_by_key(|(_, cluster)| simfunc(&cluster.centroid_seeds, kmer_freqs))
                .unwrap()
                .0;
            output[idx].push(myrseq);
        }

        output
    }
}

struct Cluster<'a> {
    centroid_seeds: SFVec,
    elements: Vec<&'a SFVec>,
}

impl<'a> Cluster<'a> {
    fn new(seeds: &'a SFVec) -> Self {
        Self { centroid_seeds: seeds.clone(), elements: vec![seeds] }
    }

    fn recenter(self) -> Self {
        if self.elements.is_empty() {
            return Self { centroid_seeds: self.centroid_seeds, elements: Vec::new() };
        }
        let new_centroid_seeds: SFVec = Clusterer::compute_cluster_centroid(&self.elements[..]);

        Self { centroid_seeds: new_centroid_seeds, elements: Vec::new() }
    }
}
