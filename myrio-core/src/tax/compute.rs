// Imports
use indicatif::MultiProgress;
use itertools::Itertools;
use rand::{SeedableRng, rngs::SmallRng, seq::IndexedRandom};

use super::Error;
use crate::{
    data::{Float, SFVec},
    tax::{
        compute_kmer_counts_for_fasta_seq, core::TaxTreeCore, kmer_store_counts_to_kmer_counts,
        store::TaxTreeStore,
    },
};

#[derive(Clone)]
pub struct TaxTreeCompute {
    pub(crate) core: TaxTreeCore<(), SFVec>,
    pub(crate) kmer_normalized_counts_fingerprint: SFVec,
    pub(crate) k_search: usize,
}

pub enum CacheOptions {
    Enabled { zstd_compression_level: i32, zstd_multithreading_opt: Option<u32> },
    Disabled,
}

impl TaxTreeCompute {
    pub fn get_kmer_normalized_counts_fingerprint(&self) -> &SFVec {
        &self.kmer_normalized_counts_fingerprint
    }

    #[allow(clippy::too_many_arguments)]
    pub fn from_store_tree(
        mut store: TaxTreeStore,
        k_cluster: usize,
        k_search: usize,
        repr_samples: usize,
        nb_bootstrap_resamples: usize,
        cache_opt: CacheOptions,
        rng: &mut impl rand::Rng,
        multi: Option<&MultiProgress>,
    ) -> Result<Self, Error> {
        let kmer_normalized_counts_fingerprint = (0..repr_samples)
            .map(|_| {
                let seq = &store.core.payloads.choose(rng).unwrap().seq;
                compute_kmer_counts_for_fasta_seq(
                    seq,
                    k_cluster,
                    nb_bootstrap_resamples,
                    &mut SmallRng::from_os_rng(),
                )
            })
            .fold(SFVec::new(0), |acc, x| acc + x)
            .into_normalized_l2();

        #[rustfmt::skip]
        let sfvec_idx = store
            .k_precomputed
            .iter()
            .copied()
            .position(|k| k_search == k)
            .map_or_else(
                || match cache_opt {
                    CacheOptions::Enabled { zstd_compression_level, zstd_multithreading_opt } => {
                        store.compute_and_append_kmer_counts(k_search, nb_bootstrap_resamples, multi);
                        store.encode_to_file(zstd_compression_level, zstd_multithreading_opt, multi)?;
                        Ok::<usize, Error>(store.k_precomputed.len())
                    },
                    CacheOptions::Disabled => {
                        store.compute_and_overwrite_kmer_counts(k_search, nb_bootstrap_resamples, multi);
                        Ok::<usize, Error>(0_usize)
                    }
                },
                Ok::<usize, Error>
        )?;

        let payloads = store
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
                gene: store.core.gene,
                highest_rank: store.core.highest_rank,
                roots: store.core.roots,
                payloads: payloads.into_boxed_slice(),
            },
            kmer_normalized_counts_fingerprint,
            k_search,
        })
    }
}
