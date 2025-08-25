use itertools::Itertools;
use rayon::{iter::ParallelIterator, slice::ParallelSlice};

use crate::{
    data::{MyrSeq, SFVec},
    similarity::{SimFunc, SimScore},
};

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

pub fn cluster(
    myrseqs: Vec<MyrSeq>,
    nb_clusters: usize,
    k: usize,
    t1: f64,
    similarity_function: SimFunc,
) -> Vec<Vec<MyrSeq>> {
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
            let new_centroid_seeds: SFVec = compute_cluster_centroid(&self.elements[..]);

            Self { centroid_seeds: new_centroid_seeds, elements: Vec::new() }
        }
    }

    // Step 1, compute the k-mer counts for each myrseq
    let (myrseqs, kcounts, nb_hcks): (Vec<MyrSeq>, Vec<SFVec>, Vec<usize>) = myrseqs
        .into_iter()
        .filter_map(|myrseq| {
            let (kcount, nb_hck) = myrseq.compute_sparse_kmer_counts(k, t1);
            if nb_hck != 0 { Some((myrseq, kcount, nb_hck)) } else { None }
        })
        .sorted_by(|(.., a), (.., b)| b.cmp(a)) // largest first
        .multiunzip();

    // Step 2, initialize the clusters
    let mut clusters: Vec<Cluster> = vec![Cluster::new(&kcounts[0])];
    let max_nb_hck = nb_hcks[0] as f64;
    while clusters.len() < nb_clusters {
        let new_cluster_seeds = kcounts
            .iter()
            .zip_eq(nb_hcks.iter())
            .min_by_key(|&(kcount, nb_hck)| {
                let mut score: f64 = clusters
                    .iter()
                    .map(|cluster| *similarity_function(&cluster.centroid_seeds, kcount))
                    .sum();

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

    // Step 3
    for _ in 0..3 {
        for kcount in kcounts.iter() {
            clusters
                .iter_mut()
                .max_by_key(|cluster| similarity_function(&cluster.centroid_seeds, kcount))
                .unwrap()
                .elements
                .push(kcount);
        }

        clusters = clusters.into_iter().map(|cluster| cluster.recenter()).collect_vec();
    }

    // Step 4
    let mut output: Vec<Vec<MyrSeq>> = vec![Vec::new(); nb_clusters];
    for (myrseq, kcount) in myrseqs.into_iter().zip_eq(kcounts.iter()) {
        let idx = clusters
            .iter()
            .enumerate()
            .max_by_key(|(_, cluster)| similarity_function(&cluster.centroid_seeds, kcount))
            .unwrap()
            .0;
        output[idx].push(myrseq);
    }

    output
}
