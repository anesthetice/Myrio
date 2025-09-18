// Imports
use indicatif::MultiProgress;
use itertools::Itertools;
use rand::{SeedableRng, seq::IndexedRandom};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use super::Error;
use crate::{
    data::{Float, SFVec},
    tax::{clade::Rank, compute_kmer_counts_for_fasta_seq, core::TaxTreeCore, store::TaxTreeStore},
};

#[derive(Clone)]
pub struct TaxTreeCompute {
    pub core: TaxTreeCore<(), SFVec>,
    /// normalized k-mer counts that represent the entire gene
    pub global_fingerprint: SFVec,
    /// set of normalized k-mer counts that represent the gene for clades at a specific rank
    pub local_fingerprints: Vec<SFVec>,
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy)]
pub enum CacheOptions {
    Enabled { zstd_compression_level: i32 },
    Disabled,
}

impl TaxTreeCompute {
    #[allow(clippy::too_many_arguments)]
    pub fn from_store_tree(
        mut ttstore: TaxTreeStore,
        cluster_k: usize,
        fingerprint_local_rank: Rank,
        fingerprint_local_nb_subsamples: usize,
        fingerprint_local_fasta_nb_resamples: usize,
        fasta_nb_resamples: usize,
        search_k: usize,
        cache_opt: CacheOptions,
        multi: Option<&MultiProgress>,
    ) -> Result<Self, Error> {
        if matches!(fingerprint_local_rank, Rank::Species) {
            return Err(Error::new_misc("Fingerprint rank cannot be set to `Species`"));
        }

        if fingerprint_local_rank > ttstore.core.highest_rank {
            return Err(Error::new_misc(format!(
                "The desired fingerprint rank of `{fingerprint_local_rank}` is higher than the tree's highest rank of `{}`",
                ttstore.core.highest_rank
            )));
        }

        let local_fingerprints: Vec<SFVec> = {
            let branches = ttstore.core.gather_branches_at_rank(fingerprint_local_rank);
            let mut output: Vec<SFVec> = Vec::with_capacity(branches.len());

            branches
                .into_par_iter()
                .map_init(rand::rngs::SmallRng::from_os_rng, |rng, branch| {
                    let leaves = branch.gather_leaves();
                    (0..fingerprint_local_nb_subsamples)
                        .map(|_| {
                            let random_payload_id = leaves.choose(rng).unwrap().payload_id;
                            let seq = &ttstore.core.payloads[random_payload_id].seq;
                            compute_kmer_counts_for_fasta_seq(
                                seq,
                                cluster_k,
                                fingerprint_local_fasta_nb_resamples,
                                rng,
                            )
                        })
                        .sum::<SFVec>()
                        .into_normalized_l2()
                })
                .collect_into_vec(&mut output);
            output
        };

        let global_fingerprint: SFVec = local_fingerprints.iter().sum::<SFVec>().into_normalized_l2();

        #[rustfmt::skip]
        let sfvec_idx = ttstore
            .k_precomputed
            .iter()
            .copied()
            .position(|k| search_k == k)
            .map_or_else(
                || match cache_opt {
                    CacheOptions::Enabled { zstd_compression_level } => {
                        ttstore.compute_and_append_kmer_counts(search_k, fasta_nb_resamples, multi);
                        ttstore.save_to_file(zstd_compression_level, multi)?;
                        Ok::<usize, Error>(ttstore.k_precomputed.len() - 1)
                    },
                    CacheOptions::Disabled => {
                        ttstore.compute_and_overwrite_kmer_counts(search_k, fasta_nb_resamples, multi);
                        Ok::<usize, Error>(0_usize)
                    }
                },
                Ok::<usize, Error>
        )?;

        let payloads = ttstore
            .core
            .payloads
            .into_iter()
            .map(|mut sp| {
                let (svec, _) = sp.kmer_store_counts_vec.swap_remove(sfvec_idx);
                svec.apply_into(Float::from).into_normalized_l2()
            })
            .collect_vec();

        Ok(Self {
            core: TaxTreeCore {
                gene: ttstore.core.gene,
                highest_rank: ttstore.core.highest_rank,
                roots: ttstore.core.roots,
                payloads: payloads.into_boxed_slice(),
            },
            global_fingerprint,
            local_fingerprints,
        })
    }
}
