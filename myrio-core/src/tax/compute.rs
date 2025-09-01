// Imports
use indicatif::MultiProgress;
use itertools::Itertools;
use rand::seq::IndexedRandom;

use super::Error;
use crate::{
    data::SFVec,
    tax::{compute_sparse_kmer_counts_for_fasta_seq, core::TaxTreeCore, store::TaxTreeStore},
};

pub struct TaxTreeCompute {
    pub(crate) core: TaxTreeCore<(), SFVec>,
    pub(crate) kmer_normcounts_repr: SFVec,
    pub(crate) k_search: usize,
}

pub enum CacheOptions {
    Enabled { zstd_compression_level: i32, zstd_multithreading_opt: Option<u32> },
    Disabled,
}

impl TaxTreeCompute {
    pub fn get_kmer_normcounts_repr(&self) -> SFVec {
        self.kmer_normcounts_repr.clone()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn from_store_tree(
        mut store: TaxTreeStore,
        k_cluster: usize,
        k_search: usize,
        repr_samples: usize,
        max_consecutive_N_before_gap: usize,
        cache_opt: CacheOptions,
        rng: &mut impl rand::Rng,
        multi: Option<&MultiProgress>,
    ) -> Result<Self, Error> {
        let mut kmer_normcounts_repr =
            if let Some(sfvec_idx) = store.pre_computed.iter().copied().position(|k| k_cluster == k) {
                (0..repr_samples)
                    .map(|_| unsafe {
                        store.core.payloads.choose(rng).unwrap().pre_comp.get_unchecked(sfvec_idx)
                    })
                    .fold(SFVec::new(0), |acc, x| acc + x)
            } else {
                (0..repr_samples)
                    .map(|_| {
                        let seq = &store.core.payloads.choose(rng).unwrap().seq;
                        compute_sparse_kmer_counts_for_fasta_seq(seq, k_cluster, max_consecutive_N_before_gap)
                    })
                    .fold(SFVec::new(0), |acc, x| acc + x)
            };
        kmer_normcounts_repr.normalize_l2();

        #[rustfmt::skip]
        let sfvec_idx = store
            .pre_computed
            .iter()
            .copied()
            .position(|k| k_search == k)
            .map_or_else(
                || match cache_opt {
                    CacheOptions::Enabled { zstd_compression_level, zstd_multithreading_opt } => {
                        store.compute_and_append_kmer_counts(k_search, max_consecutive_N_before_gap, multi);
                        store.encode_to_file(zstd_compression_level, zstd_multithreading_opt, multi)?;
                        Ok::<usize, Error>(store.pre_computed.len())
                    },
                    CacheOptions::Disabled => {
                        store.compute_and_overwrite_kmer_counts(k_search, max_consecutive_N_before_gap, multi);
                        Ok::<usize, Error>(0_usize)
                    }
                },
                Ok::<usize, Error>
        )?;

        let payloads = store
            .core
            .payloads
            .into_iter()
            .map(|mut sp| sp.pre_comp.swap_remove(sfvec_idx))
            .collect_vec();

        Ok(Self {
            core: TaxTreeCore {
                gene: store.core.gene,
                highest_rank: store.core.highest_rank,
                roots: store.core.roots,
                payloads: payloads.into_boxed_slice(),
            },
            kmer_normcounts_repr,
            k_search,
        })
    }
}
