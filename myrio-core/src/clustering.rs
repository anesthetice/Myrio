use itertools::Itertools;
use rayon::{
    iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator},
    slice::ParallelSlice,
};

use crate::{
    data::{MyrSeq, SFVec},
    similarity::{SimFunc, SimScore},
};

pub struct ClusteringParameters {
    k: usize,
    t1: f64,
    t2: usize,
    simfunc: SimFunc,
    cim: ClusterInitializationMethod,
    multithreading_flag: bool,
}

impl ClusteringParameters {
    pub fn new(
        k: usize,
        t1: f64,
        t2: usize,
        simfunc: SimFunc,
        cim: ClusterInitializationMethod,
        multithreading_flag: bool,
    ) -> Self {
        Self { k, t1, t2, simfunc, cim, multithreading_flag }
    }

    pub fn unshell(self) -> (usize, f64, usize, SimFunc, ClusterInitializationMethod, bool) {
        (self.k, self.t1, self.t2, (self.simfunc), self.cim, self.multithreading_flag)
    }
}

pub enum ClusterInitializationMethod {
    FromCentroids(Vec<SFVec>),
    FromNumber(usize),
}

pub struct ClusteringResults {
    pub clusters: Vec<Vec<MyrSeq>>,
    pub rejected: Vec<MyrSeq>,
}

pub fn compute_cluster_centroid_with_parallelism<T>(elements: &[T]) -> SFVec
where
    T: AsRef<SFVec> + Sync,
{
    let par_chunks = elements.par_chunks_exact(2);
    let remainder = par_chunks.remainder().first();

    // We have to use a slightly annoying trick as the `reduce_with` method expects the exact same input/output type and for the first 'round' input will be `&SFVec` meanwhile output will be `SFVec`
    let sfvec = par_chunks
        .map(|chunk| unsafe { chunk.get_unchecked(0).as_ref() + chunk.get_unchecked(1) })
        .reduce_with(|a, b| a + b)
        .unwrap();

    if let Some(rem) = remainder { (sfvec + rem).normalize_l2_alt() } else { sfvec.normalize_l2_alt() }
}

#[rustfmt::skip]
pub fn compute_cluster_centroid_without_parallelism<T>(elements: &[T]) -> SFVec
where
    T: AsRef<SFVec>,
{
    elements
        .iter()
        .fold(SFVec::new(0), |acc, x| acc + x)
        .normalize_l2_alt()
}

pub fn cluster(
    myrseqs: Vec<MyrSeq>,
    params: ClusteringParameters,
) -> ClusteringResults {
    let (k, t1, t2, simfunc, cim, multithreading_flag) = params.unshell();

    let mut rejected: Vec<MyrSeq> = Vec::new();

    // Step 1, compute the k-mer counts (normalized) of each myrseq
    let (myrseqs, kmer_normcounts_vec, nb_hck_vec): (Vec<MyrSeq>, Vec<SFVec>, Vec<usize>) = myrseqs
        .into_iter()
        .filter_map(|myrseq| {
            let (kmer_normcounts, nb_hck) = myrseq.compute_sparse_kmer_normcounts(k, t1);
            if nb_hck > t2 {
                Some((myrseq, kmer_normcounts, nb_hck))
            } else {
                rejected.push(myrseq);
                None
            }
        })
        //.sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
        .multiunzip();

    // Step 2, initialize the clusters
    let mut clusters = match cim {
        ClusterInitializationMethod::FromCentroids(centroids) => centroids
            .into_iter()
            .map(|centroid| Cluster { centroid_seeds: centroid, elements: Vec::new() })
            .collect_vec(),
        ClusterInitializationMethod::FromNumber(nb_clusters) => {
            let mut clusters: Vec<Cluster> = vec![Cluster::new(&kmer_normcounts_vec[0])];
            let max_nb_hck = nb_hck_vec[0] as f64;
            while clusters.len() < nb_clusters {
                let new_cluster_seeds = kmer_normcounts_vec
                    .iter()
                    .zip_eq(nb_hck_vec.iter())
                    .min_by_key(|&(kcount, nb_hck)| {
                        let mut score: f64 =
                            clusters.iter().map(|cluster| *simfunc(&cluster.centroid_seeds, kcount)).sum();
                        // penalty (increases score) for seqs with a low number of high-confidence k-mers
                        score = 0.7 * score + 0.3 * (1.0 - (*nb_hck as f64) / max_nb_hck);
                        SimScore::try_new(score).unwrap()
                    })
                    .unwrap()
                    .0;
                clusters.push(Cluster::new(new_cluster_seeds));
            }
            clusters
        }
    };

    // Step 3
    if multithreading_flag {
        let mut target: Vec<(&SFVec, usize)> = Vec::with_capacity(kmer_normcounts_vec.len());
        for _ in 0..3 {
            kmer_normcounts_vec
                .par_iter()
                .map(|kmer_normcounts| {
                    let cl_idx = clusters
                        .iter()
                        .enumerate()
                        .max_by_key(|(_, cluster)| simfunc(&cluster.centroid_seeds, kmer_normcounts))
                        .unwrap()
                        .0;
                    (kmer_normcounts, cl_idx)
                })
                .collect_into_vec(&mut target);
            target.iter().for_each(|(sfvec, cl_idx)| clusters[*cl_idx].elements.push(sfvec));
            target.clear();
            clusters = clusters.into_iter().map(Cluster::recenter_with_parallelism).collect_vec();
        }
    } else {
        for _ in 0..3 {
            kmer_normcounts_vec.iter().for_each(|kmer_normcounts| {
                clusters
                    .iter_mut()
                    .max_by_key(|cluster| simfunc(&cluster.centroid_seeds, kmer_normcounts))
                    .unwrap()
                    .elements
                    .push(kmer_normcounts);
            });
            clusters = clusters.into_iter().map(Cluster::recenter_without_parallelism).collect_vec();
        }
    }

    // Step 4
    // (could also add multi-threading to this part in the future)
    let mut output: Vec<Vec<MyrSeq>> = vec![Vec::new(); clusters.len()];
    for (myrseq, kmer_normcounts) in myrseqs.into_iter().zip_eq(kmer_normcounts_vec.iter()) {
        let idx = clusters
            .iter()
            .enumerate()
            .max_by_key(|(_, cluster)| simfunc(&cluster.centroid_seeds, kmer_normcounts))
            .unwrap()
            .0;
        output[idx].push(myrseq);
    }

    ClusteringResults { clusters: output, rejected }
}

struct Cluster<'a> {
    centroid_seeds: SFVec,
    elements: Vec<&'a SFVec>,
}

impl<'a> Cluster<'a> {
    fn new(seeds: &'a SFVec) -> Self {
        Self { centroid_seeds: seeds.clone(), elements: vec![seeds] }
    }

    fn recenter_with_parallelism(self) -> Self {
        if self.elements.is_empty() {
            return Self { centroid_seeds: self.centroid_seeds, elements: Vec::new() };
        }
        let new_centroid_seeds: SFVec = compute_cluster_centroid_with_parallelism(&self.elements[..]);

        Self { centroid_seeds: new_centroid_seeds, elements: Vec::new() }
    }

    fn recenter_without_parallelism(self) -> Self {
        if self.elements.is_empty() {
            return Self { centroid_seeds: self.centroid_seeds, elements: Vec::new() };
        }
        let new_centroid_seeds: SFVec = compute_cluster_centroid_without_parallelism(&self.elements[..]);

        Self { centroid_seeds: new_centroid_seeds, elements: Vec::new() }
    }
}
