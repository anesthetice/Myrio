#![allow(non_snake_case)]

use rayon::iter::IntoParallelIterator;

// Imports
use crate::{
    data::SFVec,
    similarity::SimScore,
    tax::{Error, clade::Rank, core::TaxTreeCore, store::TaxTreeStore},
};

pub struct TaxTreeResults {
    core: TaxTreeCore<(), SimScore>,
}

impl TaxTreeResults {
    pub fn from_store_with_caching(
        mut store: TaxTreeStore,
        k: usize,
        max_consecutive_N_before_gap: usize,
        compression_level: i32,
        multithreading_flag: bool,
        nb_threads_available: usize,
    ) -> Result<Self, Error> {
        #[rustfmt::skip]
        let sfvec_idx = store
            .pre_computed
            .iter()
            .copied()
            .find(|&k_| k == k_)
            .map_or_else(
                || {
                    store.compute_and_append_kmer_counts(k, max_consecutive_N_before_gap);
                    store.encode_to_file(compression_level, multithreading_flag, nb_threads_available)?;
                    Ok::<usize, Error>(store.pre_computed.len())
                },
                Ok::<usize, Error>
        )?;

        Ok(Self::from_store_core(store, k, sfvec_idx))
    }

    pub fn from_store_without_caching(
        mut store: TaxTreeStore,
        k: usize,
        max_consecutive_N_before_gap: usize,
    ) -> Self {
        #[rustfmt::skip]
        let sfvec_idx = store
            .pre_computed
            .iter()
            .copied()
            .find(|&k_| k == k_)
            .unwrap_or_else(|| {
                store.compute_and_overwrite_kmer_counts(k, max_consecutive_N_before_gap);
                0_usize
            });
        Self::from_store_core(store, k, sfvec_idx)
    }

    fn from_store_core(
        mut store: TaxTreeStore,
        k: usize,
        sfvec_idx: usize,
    ) -> Self {
        /*
        fn dive_recursive(
            store_node: StoreNode,
            above: &mut Vec<Node>,
            sfvec_idx: usize,
        ) {
            match store_node {
                StoreNode::Branch { name, children } => {
                    let mut current: Vec<Node> = Vec::with_capacity(children.len());
                    for store_node in children {
                        dive_recursive(store_node, &mut current, sfvec_idx);
                    }
                    above.push(Node::new_branch(name, current.into_boxed_slice()));
                }
                StoreNode::Leaf { name, seq: _, mut pre_comp } => {
                    above.push(Node::new_leaf(name, pre_comp.swap_remove(sfvec_idx)))
                }
            }
        }
        */
        store.core.payloads.into_par_iter();
        todo!()
    }
}
