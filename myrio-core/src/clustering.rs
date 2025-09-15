// Imports
use std::ops::Div;

use console::style;
use itertools::Itertools;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};

use crate::{
    data::{Float, MyrSeq, SFVec},
    similarity::{SimFunc, SimScore, Similarity},
};

#[derive(Debug, Clone)]
pub struct ClusteringParameters {
    /// k-mer length, `6` seems to be the most sensible value for clustering
    pub k: usize,
    /// Score cutoff threshold i.e., if the product of error probabilities of the nucleotides in a k-mer is below this threshold, then it is ignored
    pub t1: Float,
    /// Minimum number of valid k-mers for a sequence to not be discarded
    pub t2: usize,
    /// The type of similarity to use
    pub similarity: Similarity,
    /// Can be empty, should be normalized
    pub intial_centroids: Vec<SFVec>,
    /// Total number of clusters expected, should logically be greater or equal to `initial_centroids.len()`
    pub expected_nb_of_clusters: usize,
    /// Fraction to stop the main clustering iteration (not recommended to be above `1E-2`)
    pub eta_improvement: Float,
    /// Maximum number of iterations before stopping
    pub nb_iters_max: usize,
    /// If set to `Some(ssdcf_factor)`, then within the semi-final clusters, any sequence for which the silhouette score is below `mean - ssdcf_factor * std` will be discarded
    /// Note that `ssdcf_factor` stands for "silhouette std deviation cutoff factor"
    pub silhouette_trimming: Option<Float>,
}

pub struct ClusteringOutput {
    pub clusters: Vec<Vec<MyrSeq>>,
    pub rejected: Vec<MyrSeq>,
}

/// The input `members` should ideally not be normalized, output is normalized
pub fn compute_cluster_centroid<T>(members: &[T]) -> SFVec
where
    T: AsRef<SFVec> + Sync,
{
    members
        .par_iter()
        .fold(|| SFVec::new(0), |a, b| a + b.as_ref())
        .reduce(|| SFVec::new(0), |a, b| a + b)
        .into_normalized_l2()
}

pub fn cluster(
    myrseqs: Vec<MyrSeq>,
    params: ClusteringParameters,
) -> ClusteringOutput {
    let ClusteringParameters {
        k,
        t1,
        t2,
        similarity,
        intial_centroids,
        expected_nb_of_clusters,
        eta_improvement,
        nb_iters_max,
        silhouette_trimming,
    } = params;

    if intial_centroids.is_empty() && expected_nb_of_clusters == 0 {
        panic!("Ensure `initial_centroids` is not empty or `expected_nb_of_clusters` > 0")
    }

    let mut rejected: Vec<MyrSeq> = Vec::new();

    // Step 1, compute the k-mer counts of each myrseq, discarding k-mers with low quality then discarding sequences with not enough remaining valid k-mers
    let (myrseqs, kmer_counts_vec, nb_hck_vec): (Vec<MyrSeq>, Vec<SFVec>, Vec<usize>) = myrseqs
        .into_iter()
        .filter_map(|myrseq| {
            let (kmer_counts, nb_hck) = myrseq.compute_kmer_counts(k, t1);
            if nb_hck > t2 {
                Some((myrseq, kmer_counts, nb_hck))
            } else {
                rejected.push(myrseq);
                None
            }
        })
        //.sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
        .multiunzip();

    // We keep both the non-normalized k-mer counts (for cluster members) but also the normalized ones (for more efficient computation)
    let kmer_counts_normalized_vec =
        kmer_counts_vec.clone().into_iter().map(SFVec::into_normalized_l2).collect_vec();

    let simfunc: SimFunc = similarity.to_simfunc(true);

    // Step 2, initialize the clusters
    let mut clusters = intial_centroids
        .into_iter()
        .map(|centroid| Cluster { centroid, members: Vec::new(), members_normalized: Vec::new() })
        .collect_vec();

    let (max_nb_hck_idx, max_nb_hck) =
        nb_hck_vec.iter().enumerate().max_by_key(|(_, nb_hck)| **nb_hck).unwrap();
    let max_nb_hck = *max_nb_hck as Float;

    if clusters.is_empty() {
        clusters.push(Cluster::from_normalized_centroid(&kmer_counts_normalized_vec[max_nb_hck_idx]));
    }

    while clusters.len() < expected_nb_of_clusters {
        let centroid_for_new_cluster = kmer_counts_normalized_vec
            .iter()
            .zip_eq(nb_hck_vec.iter())
            .min_by_key(|&(kmer_counts_normalized, nb_hck)| {
                let mut score = clusters
                    .iter()
                    .map(|cluster| *simfunc(&cluster.centroid, kmer_counts_normalized))
                    .sum::<Float>()
                    .div(clusters.len() as Float);
                // add a penalty (increases score) for sequences with a low number of high-confidence k-mers
                score *= 0.8 + 0.2 * (1.0 - (*nb_hck as Float) / max_nb_hck);
                SimScore::try_new(score).unwrap()
            })
            .unwrap()
            .0;
        clusters.push(Cluster::from_normalized_centroid(centroid_for_new_cluster));
    }

    // Step 3, main clustering loop
    let mut target: Vec<(usize, usize)> = Vec::with_capacity(kmer_counts_normalized_vec.len());
    let mut previous_mean_wcsss: Float = 1E-6;
    let mut improvement_score: Float = Float::MAX;
    let mut nb_iters: usize = 0;

    while improvement_score > eta_improvement && nb_iters < nb_iters_max {
        kmer_counts_normalized_vec
            .par_iter()
            .enumerate()
            .map(|(kc_idx, kmer_counts_normalized)| {
                let cl_idx = clusters
                    .iter()
                    .enumerate()
                    .max_by_key(|(_, cluster)| simfunc(&cluster.centroid, kmer_counts_normalized))
                    .unwrap()
                    .0;
                (kc_idx, cl_idx)
            })
            .collect_into_vec(&mut target);

        target.iter().for_each(|&(kc_idx, cl_idx)| unsafe {
            clusters
                .get_unchecked_mut(cl_idx)
                .push(kmer_counts_vec.get_unchecked(kc_idx), kmer_counts_normalized_vec.get_unchecked(kc_idx))
        });
        target.clear();

        let current_mean_wcsss =
            clusters.iter().map(|cl| cl.wcsss(simfunc)).sum::<Float>() / expected_nb_of_clusters as Float;
        improvement_score = (current_mean_wcsss - previous_mean_wcsss) / previous_mean_wcsss;
        previous_mean_wcsss = current_mean_wcsss;

        #[cfg(debug_assertions)]
        println!("improvement score: {improvement_score}");

        clusters = clusters.into_iter().map(Cluster::recenter).collect_vec();
        nb_iters += 1;
    }

    /*
    {
        kmer_counts_normalized_vec
            .par_iter()
            .enumerate()
            .map(|(kc_idx, kmer_counts_normalized)| {
                let cl_idx = clusters
                    .iter()
                    .enumerate()
                    .max_by_key(|(_, cluster)| simfunc(&cluster.centroid, kmer_counts_normalized))
                    .unwrap()
                    .0;
                (kc_idx, cl_idx)
            })
            .collect_into_vec(&mut target);

        target.iter().for_each(|&(kc_idx, cl_idx)| unsafe {
            clusters
                .get_unchecked_mut(cl_idx)
                .push(kmer_counts_vec.get_unchecked(kc_idx), kmer_counts_normalized_vec.get_unchecked(kc_idx))
        });
        target.clear();

        Cluster::print_analyze(&clusters, None, similarity, &["ITS", "matK", "rbcL"]);
        clusters = clusters.into_iter().map(Cluster::recenter).collect_vec();
    }
    */

    // Step 4 bypass if silhouette_trimming is set to `None` (i.e., disabled)
    if silhouette_trimming.is_none() {
        let mut output: Vec<Vec<MyrSeq>> = vec![Vec::new(); clusters.len()];
        for (myrseq, kmer_counts_normalized) in myrseqs.into_iter().zip_eq(kmer_counts_normalized_vec.iter())
        {
            let idx = clusters
                .iter()
                .enumerate()
                .max_by_key(|(_, cluster)| simfunc(&cluster.centroid, kmer_counts_normalized))
                .unwrap()
                .0;
            output[idx].push(myrseq);
        }
        return ClusteringOutput { clusters: output, rejected };
    }

    // "silhouette std deviation cutoff factor"
    let ssdcf = silhouette_trimming.unwrap();

    // Step 4 and 5, silhouette thinning and simply computing the output
    // (could also add multi-threading to this part in the future)
    let mut clustered_myrseqs: Vec<Vec<MyrSeq>> = vec![Vec::new(); clusters.len()];
    for ((myrseq, kmer_counts_normalized), kmer_counts) in myrseqs
        .into_iter()
        .zip_eq(kmer_counts_normalized_vec.iter())
        .zip_eq(kmer_counts_vec.iter())
    {
        let (cl_idx, cluster) = clusters
            .iter_mut()
            .enumerate()
            .max_by_key(|(_, cluster)| simfunc(&cluster.centroid, kmer_counts_normalized))
            .unwrap();

        cluster.push(kmer_counts, kmer_counts_normalized);
        unsafe { clustered_myrseqs.get_unchecked_mut(cl_idx).push(myrseq) };
    }

    let silhouette_scores_vec = Cluster::compute_silhouette_scores_vec(&clusters, similarity);

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

            let left_cutoff = mean - std * ssdcf;

            #[cfg(debug_assertions)]
            println!("mean={mean}, std={std}, left_cutoff={left_cutoff:.3}");

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

    ClusteringOutput { clusters: output, rejected }
}

struct Cluster<'a> {
    /// Should be normalized
    centroid: SFVec,
    members: Vec<&'a SFVec>,
    members_normalized: Vec<&'a SFVec>,
}

impl<'a> Cluster<'a> {
    fn from_normalized_centroid(centroid: &'a SFVec) -> Self {
        Self {
            centroid: centroid.clone().into_normalized_l2(),
            members: Vec::new(),
            members_normalized: Vec::new(),
        }
    }

    fn recenter(self) -> Self {
        if self.members.is_empty() {
            return Self { centroid: self.centroid, members: Vec::new(), members_normalized: Vec::new() };
        }
        let new_centroid: SFVec = compute_cluster_centroid(&self.members[..]);

        Self { centroid: new_centroid, members: Vec::new(), members_normalized: Vec::new() }
    }

    fn push(
        &mut self,
        kmer_counts: &'a SFVec,
        kmer_counts_normalized: &'a SFVec,
    ) {
        self.members.push(kmer_counts);
        self.members_normalized.push(kmer_counts_normalized);
    }

    /// Within-cluster sum of similarity scores
    fn wcsss(
        &self,
        simfunc: SimFunc,
    ) -> Float {
        self.members_normalized
            .iter()
            .map(|&sfvec| *simfunc(sfvec, &self.centroid))
            .sum1::<Float>()
            .unwrap_or_default()
    }

    fn compute_silhouette_scores_vec(
        clusters: &[Self],
        similarity: Similarity,
    ) -> Vec<Vec<Float>> {
        #[rustfmt::skip]
        fn inner(cluster: &Cluster, neighbor: &Cluster, index: usize, simfunc: SimFunc) -> Float {
            let i = cluster.members_normalized[index];

            let a = cluster.members_normalized.iter()
                .map(|&e| *Similarity::dist(simfunc(i, e)))
                .sum::<Float>()
                .div((cluster.members_normalized.len() - 1) as Float);
            let b = neighbor.members_normalized.iter()
                .map(|&e| *Similarity::dist(simfunc(i, e)))
                .sum::<Float>()
                .div(neighbor.members_normalized.len() as Float);

            (b - a) / a.max(b)
        }

        let simfunc: SimFunc = similarity.to_simfunc(true);
        let mut output: Vec<Vec<Float>> = Vec::with_capacity(clusters.len());

        clusters
            .par_iter()
            .enumerate()
            .map(|(i, cluster)| {
                let neighbor = clusters
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .max_by_key(|(_, other_cluster)| simfunc(&cluster.centroid, &other_cluster.centroid))
                    .unwrap()
                    .1;

                let mut silscores: Vec<Float> = vec![0.0; cluster.members.len()];

                silscores
                    .par_iter_mut()
                    .enumerate()
                    .for_each(|(index, value)| *value = inner(cluster, neighbor, index, simfunc));

                silscores
            })
            .collect_into_vec(&mut output);

        output
    }

    #[allow(unused)]
    fn print_analyze(
        clusters: &'a [Cluster],
        silhouette_scores_vec: Option<&[Vec<Float>]>,
        similarity: Similarity,
        names: &[&str],
    ) {
        let silhouette_scores_vec_computed = if silhouette_scores_vec.is_none() {
            Self::compute_silhouette_scores_vec(clusters, similarity)
        } else {
            Vec::with_capacity(0)
        };

        let silhouette_scores_vec = if let Some(silscores_vec) = silhouette_scores_vec {
            silscores_vec
        } else {
            silhouette_scores_vec_computed.as_slice()
        };

        let simfunc: SimFunc = similarity.to_simfunc(true);

        debug_assert_eq!(clusters.len(), silhouette_scores_vec.len(),);

        debug_assert_eq!(clusters.len(), names.len());

        for (cluster, silscores, name) in itertools::izip!(clusters, silhouette_scores_vec, names) {
            let n = silscores.len() as Float;

            let (wcsss_weighted, mean, std) = if !silscores.is_empty() {
                let wcss_weighted: Float = cluster.wcsss(simfunc) / n;
                let mean: Float = silscores.iter().sum::<Float>() / n;
                let std: Float =
                    (n.powi(-1) * silscores.iter().map(|x| (x - mean).powi(2)).sum::<Float>()).sqrt();
                (wcss_weighted, mean, std)
            } else {
                (0.0, 0.0, 0.0)
            };
            let name_bold = style(name).bold();
            indoc::printdoc! {"
                {name_bold} cluster
                ├── within-cluster sum of similarity scores weighted: {wcsss_weighted:.4}
                └── silhouette: (μ={mean:.3}, σ={std:.3})
            "};
        }
    }
}
