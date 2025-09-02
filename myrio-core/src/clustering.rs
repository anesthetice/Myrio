use std::{
    f64,
    ops::{Div, Sub},
};

use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
    data::{Float, MyrSeq, SFVec},
    similarity::{SimFunc, SimScore, Similarity},
};

#[derive(Debug, Clone)]
pub struct ClusteringParameters {
    pub k: usize,
    pub t1: f64,
    pub t2: usize,
    pub similarity: Similarity,
    pub cim: ClusterInitializationMethod,
    pub eta_improvement: Float,
    pub nb_iters_max: usize,
    pub silhouette_std_deviation_cutoff_factor: Float,
    pub multithreading_flag: bool,
}

impl ClusteringParameters {
    pub fn new(
        k: usize,
        t1: f64,
        t2: usize,
        similarity: Similarity,
        cim: ClusterInitializationMethod,
        eta_improvement: Float,
        nb_iters_max: usize,
        silhouette_std_deviation_cutoff_factor: Float,
        multithreading_flag: bool,
    ) -> Self {
        Self {
            k,
            t1,
            t2,
            similarity,
            cim,
            eta_improvement,
            nb_iters_max,
            silhouette_std_deviation_cutoff_factor,
            multithreading_flag,
        }
    }

    pub fn unshell(
        self
    ) -> (usize, f64, usize, Similarity, ClusterInitializationMethod, Float, usize, Float, bool) {
        (
            self.k,
            self.t1,
            self.t2,
            self.similarity,
            self.cim,
            self.eta_improvement,
            self.nb_iters_max,
            self.silhouette_std_deviation_cutoff_factor,
            self.multithreading_flag,
        )
    }
}

#[derive(Debug, Clone)]
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
    elements
        .par_iter()
        .fold(|| SFVec::new(0), |a, b| a + b.as_ref())
        .reduce(|| SFVec::new(0), |a, b| a + b)
        .normalize_l2_alt()
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
    let (
        k,
        t1,
        t2,
        similarity,
        cim,
        eta_improvement,
        nb_iters_max,
        silhouette_std_deviation_cutoff_factor,
        multithreading_flag,
    ) = params.unshell();

    let mut rejected: Vec<MyrSeq> = Vec::new();

    // Step 1, compute the k-mer counts of each myrseq, normalized either way (too bothersome later on otherwise)
    let (myrseqs, kmer_norm_counts_vec, nb_hck_vec): (Vec<MyrSeq>, Vec<SFVec>, Vec<usize>) = myrseqs
        .into_iter()
        .filter_map(|myrseq| {
            let (kmer_norm_counts, nb_hck) = myrseq.compute_sparse_kmer_counts_normalized(k, t1);
            if nb_hck > t2 {
                Some((myrseq, kmer_norm_counts, nb_hck))
            } else {
                rejected.push(myrseq);
                None
            }
        })
        //.sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
        .multiunzip();
    let simfunc: SimFunc = similarity.to_simfunc(true);

    // Step 2, initialize the clusters
    let mut clusters = match cim {
        ClusterInitializationMethod::FromCentroids(centroids) => centroids
            .into_iter()
            .map(|centroid| Cluster { centroid_seeds: centroid, elements: Vec::new() })
            .collect_vec(),
        ClusterInitializationMethod::FromNumber(nb_clusters) => {
            let mut clusters: Vec<Cluster> = vec![Cluster::new(&kmer_norm_counts_vec[0])];
            let max_nb_hck = nb_hck_vec[0] as Float;
            while clusters.len() < nb_clusters {
                let new_cluster_seeds = kmer_norm_counts_vec
                    .iter()
                    .zip_eq(nb_hck_vec.iter())
                    .min_by_key(|&(kcount, nb_hck)| {
                        let mut score = clusters
                            .iter()
                            .map(|cluster| *simfunc(&cluster.centroid_seeds, kcount))
                            .sum::<Float>();
                        // penalty (increases score) for seqs with a low number of high-confidence k-mers
                        score = 0.7 * score + 0.3 * (1.0 - (*nb_hck as Float) / max_nb_hck);
                        SimScore::try_new(score).unwrap()
                    })
                    .unwrap()
                    .0;
                clusters.push(Cluster::new(new_cluster_seeds));
            }
            clusters
        }
    };
    let nb_clusters = clusters.len();

    // Step 3, main clustering loop
    if multithreading_flag {
        let mut target: Vec<(&SFVec, usize)> = Vec::with_capacity(kmer_norm_counts_vec.len());
        let mut previous_mean_wcsss: Float = 1E-6;
        let mut improvement_score: Float = Float::MAX;
        let mut nb_iters: usize = 0;

        while improvement_score > eta_improvement && nb_iters < nb_iters_max {
            kmer_norm_counts_vec
                .par_iter()
                .map(|kmer_norm_counts| {
                    let cl_idx = clusters
                        .iter()
                        .enumerate()
                        .max_by_key(|(_, cluster)| simfunc(&cluster.centroid_seeds, kmer_norm_counts))
                        .unwrap()
                        .0;
                    (kmer_norm_counts, cl_idx)
                })
                .collect_into_vec(&mut target);
            target.iter().for_each(|(sfvec, cl_idx)| clusters[*cl_idx].elements.push(sfvec));
            target.clear();

            let current_mean_wcsss =
                clusters.iter().map(|cl| cl.wcsss(simfunc)).sum::<Float>() / nb_clusters as Float;
            improvement_score = (current_mean_wcsss - previous_mean_wcsss) / previous_mean_wcsss;
            previous_mean_wcsss = current_mean_wcsss;

            clusters = clusters.into_iter().map(Cluster::recenter_without_parallelism).collect_vec();
            nb_iters += 1;
        }
    } else {
        let mut previous_mean_wcsss: Float = 1E-6;
        let mut improvement_score: Float = Float::MAX;
        let mut nb_iters: usize = 0;

        while improvement_score > eta_improvement && nb_iters < nb_iters_max {
            kmer_norm_counts_vec.iter().for_each(|kmer_norm_counts| {
                clusters
                    .iter_mut()
                    .max_by_key(|cluster| simfunc(&cluster.centroid_seeds, kmer_norm_counts))
                    .unwrap()
                    .elements
                    .push(kmer_norm_counts);
            });

            let current_mean_wcsss =
                clusters.iter().map(|cl| cl.wcsss(simfunc)).sum::<Float>() / nb_clusters as Float;
            improvement_score = (current_mean_wcsss - previous_mean_wcsss) / previous_mean_wcsss;
            previous_mean_wcsss = current_mean_wcsss;

            #[cfg(debug_assertions)]
            eprintln!("improvement score: {improvement_score}");

            clusters = clusters.into_iter().map(Cluster::recenter_without_parallelism).collect_vec();
            nb_iters += 1;
        }
    }

    // Step 4 and 5, silhouette thinning and simply computing the output
    // (could also add multi-threading to this part in the future)
    let mut clustered_myrseqs: Vec<Vec<MyrSeq>> = vec![Vec::new(); clusters.len()];
    for (myrseq, kmer_norm_counts) in myrseqs.into_iter().zip_eq(kmer_norm_counts_vec.iter()) {
        let (idx, cluster) = clusters
            .iter_mut()
            .enumerate()
            .max_by_key(|(_, cluster)| simfunc(&cluster.centroid_seeds, kmer_norm_counts))
            .unwrap();
        cluster.elements.push(kmer_norm_counts);
        clustered_myrseqs[idx].push(myrseq);
    }
    let silhouette_scores_vec = Cluster::compute_silhouette_scores_vec(&clusters, simfunc);

    let output = clustered_myrseqs
        .into_iter()
        .zip_eq(silhouette_scores_vec)
        .map(|(myrseq_cluster, silhouette_scores)| {
            if silhouette_scores.is_empty() {
                return myrseq_cluster;
            }
            let n = silhouette_scores.len() as Float;
            let mean = silhouette_scores.iter().sum::<Float>() / n;
            let std =
                (n.powi(-1) * silhouette_scores.iter().map(|x| (x - mean).powi(2)).sum::<Float>()).sqrt();

            let left_cutoff = mean - std * silhouette_std_deviation_cutoff_factor;

            myrseq_cluster
                .into_iter()
                .zip_eq(silhouette_scores)
                .filter_map(|(myrseq, silscore)| {
                    if silscore < left_cutoff {
                        rejected.push(myrseq);
                        None
                    } else {
                        Some(myrseq)
                    }
                })
                .collect_vec()
        })
        .collect_vec();

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

    /// Within-cluster sum of similarity scores
    fn wcsss(
        &self,
        simfunc: SimFunc,
    ) -> Float {
        self.elements
            .iter()
            .map(|&sfvec| *simfunc(sfvec, &self.centroid_seeds))
            .sum1::<Float>()
            .unwrap_or_default()
    }

    fn compute_silhouette_scores_vec(
        clusters: &[Self],
        simfunc: SimFunc,
    ) -> Vec<Vec<Float>> {
        #[rustfmt::skip]
        fn inner(cluster: &Cluster, neighbor: &Cluster, index: usize, simfunc: SimFunc) -> Float {
            let i = cluster.elements[index];

            let a = cluster.elements.iter()
                .map(|&e| *Similarity::dist(simfunc(i, e)))
                .sum::<Float>()
                .sub(1.0)
                .div((cluster.elements.len() - 1) as Float);
            let b = neighbor.elements.iter()
                .map(|&e| *Similarity::dist(simfunc(i, e)))
                .sum::<Float>()
                .div(neighbor.elements.len() as Float);

            (b - a) / a.max(b)
        }

        let mut output = Vec::with_capacity(clusters.len());

        for (i, cluster) in clusters.iter().enumerate() {
            let neighbor = clusters
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .max_by_key(|(_, other_cluster)| {
                    simfunc(&cluster.centroid_seeds, &other_cluster.centroid_seeds)
                })
                .unwrap()
                .1;

            output.push(
                (0..cluster.elements.len())
                    .map(|index| inner(cluster, neighbor, index, simfunc))
                    .collect_vec(),
            );
        }
        output
    }
}
