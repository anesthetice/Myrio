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
}

#[derive(Debug, Clone, Copy)]
pub enum CacheOptions {
    Enabled { zstd_compression_level: i32 },
    Disabled,
}

impl TaxTreeCompute {
    pub fn get_kmer_normalized_counts_fingerprint(&self) -> &SFVec {
        &self.kmer_normalized_counts_fingerprint
    }

    #[allow(clippy::too_many_arguments)]
    pub fn from_store_tree(
        mut ttstore: TaxTreeStore,
        cluster_k: usize,
        fingerprint_nb_subsamples: usize,
        fingerprint_fasta_nb_resamples: usize,
        search_k: usize,
        search_fasta_nb_resamples: usize,
        cache_opt: CacheOptions,
        rng: &mut impl rand::Rng,
        multi: Option<&MultiProgress>,
    ) -> Result<Self, Error> {
        let kmer_normalized_counts_fingerprint = (0..fingerprint_nb_subsamples)
            .map(|_| {
                let seq = &ttstore.core.payloads.choose(rng).unwrap().seq;
                compute_kmer_counts_for_fasta_seq(
                    seq,
                    cluster_k,
                    fingerprint_fasta_nb_resamples,
                    &mut SmallRng::from_os_rng(),
                )
            })
            .fold(SFVec::new(0), |acc, x| acc + x)
            .into_normalized_l2();

        #[rustfmt::skip]
        let sfvec_idx = ttstore
            .k_precomputed
            .iter()
            .copied()
            .position(|k| search_k == k)
            .map_or_else(
                || match cache_opt {
                    CacheOptions::Enabled { zstd_compression_level } => {
                        ttstore.compute_and_append_kmer_counts(search_k, search_fasta_nb_resamples, multi);
                        ttstore.save_to_file(zstd_compression_level, multi)?;
                        Ok::<usize, Error>(ttstore.k_precomputed.len())
                    },
                    CacheOptions::Disabled => {
                        ttstore.compute_and_overwrite_kmer_counts(search_k, search_fasta_nb_resamples, multi);
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
            kmer_normalized_counts_fingerprint,
        })
    }
}
