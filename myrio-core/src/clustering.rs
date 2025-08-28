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
    multi_threading_flag: bool,
}

impl ClusteringParameters {
    pub fn new(
        k: usize,
        t1: f64,
        t2: usize,
        simfunc: SimFunc,
        cim: ClusterInitializationMethod,
        multi_threading_flag: bool,
    ) -> Self {
        Self { k, t1, t2, simfunc, cim, multi_threading_flag }
    }

    pub fn unshell(self) -> (usize, f64, usize, SimFunc, ClusterInitializationMethod, bool) {
        (self.k, self.t1, self.t2, (self.simfunc), self.cim, self.multi_threading_flag)
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
    let n = elements.len() as f64;

    let par_chunks = elements.par_chunks_exact(2);
    let remainder = par_chunks.remainder().first();

    // We have to use a slightly annoying trick as the `reduce_with` method expects the exact same input/output type and for the first 'round' input will be `&SFVec` meanwhile output will be `SFVec`
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

pub fn compute_cluster_centroid_without_parallelism<T>(elements: &[T]) -> SFVec
where
    T: AsRef<SFVec>,
{
    let n = elements.len() as f64;

    let mut sfvec = elements
        .iter()
        .fold(SFVec::new(0), |acc, x| acc.merge_with_and_apply(x.as_ref(), |x, y| x + y));
    sfvec /= n;
    sfvec
}

pub fn cluster(
    myrseqs: Vec<MyrSeq>,
    params: ClusteringParameters,
) -> ClusteringResults {
    let (k, t1, t2, simfunc, cim, multi_threading_flag) = params.unshell();

    let mut rejected: Vec<MyrSeq> = Vec::new();

    // Step 1, compute the k-mer freqs of each myrseq
    let (myrseqs, kmer_freqs_vec, nb_hck_vec): (Vec<MyrSeq>, Vec<SFVec>, Vec<usize>) = myrseqs
        .into_iter()
        .filter_map(|myrseq| {
            let (mut kmer_freqs, nb_hck) = myrseq.compute_sparse_kmer_counts(k, t1);
            if nb_hck > t2 {
                kmer_freqs /= nb_hck as f64;
                Some((myrseq, kmer_freqs, nb_hck))
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
            let mut clusters: Vec<Cluster> = vec![Cluster::new(&kmer_freqs_vec[0])];
            let max_nb_hck = nb_hck_vec[0] as f64;
            while clusters.len() < nb_clusters {
                let new_cluster_seeds = kmer_freqs_vec
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
    if multi_threading_flag {
        let mut target: Vec<(&SFVec, usize)> = Vec::with_capacity(kmer_freqs_vec.len());
        for _ in 0..3 {
            kmer_freqs_vec
                .par_iter()
                .map(|kmer_freqs| {
                    let cl_idx = clusters
                        .iter()
                        .enumerate()
                        .max_by_key(|(_, cluster)| simfunc(&cluster.centroid_seeds, kmer_freqs))
                        .unwrap()
                        .0;
                    (kmer_freqs, cl_idx)
                })
                .collect_into_vec(&mut target);
            target.iter().for_each(|(sfvec, cl_idx)| clusters[*cl_idx].elements.push(sfvec));
            target.clear();
            clusters = clusters.into_iter().map(Cluster::recenter_with_parallelism).collect_vec();
        }
    } else {
        for _ in 0..3 {
            kmer_freqs_vec.iter().for_each(|kmer_freqs| {
                clusters
                    .iter_mut()
                    .max_by_key(|cluster| simfunc(&cluster.centroid_seeds, kmer_freqs))
                    .unwrap()
                    .elements
                    .push(kmer_freqs);
            });
            clusters = clusters.into_iter().map(Cluster::recenter_without_parallelism).collect_vec();
        }
    }

    // Step 4
    // (could also add multi-threading to this part in the future)
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
