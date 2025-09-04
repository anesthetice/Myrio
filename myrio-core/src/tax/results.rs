#![allow(non_snake_case)]

use std::{ops::Neg, usize};

use indicatif::MultiProgress;
// Imports
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use crate::data::SparseVec;
use crate::{
    data::{Float, SFVec},
    similarity::{SimFunc, SimScore, Similarity},
    tax::{
        Error,
        compute::TaxTreeCompute,
        core::{Leaf, Node, TaxTreeCore},
    },
};

#[derive(Clone, Copy)]
pub struct BranchExtra {
    pub mean: Float,
}

impl std::fmt::Display for BranchExtra {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "μ={:.3}", self.mean)
    }
}

pub struct TaxTreeResults {
    pub core: TaxTreeCore<BranchExtra, SimScore>,
}

impl TaxTreeResults {
    pub fn from_compute_tree(
        query: SFVec,
        compute: TaxTreeCompute,
        similarity: Similarity,
        multi: Option<&MultiProgress>,
    ) -> Result<Self, Error> {
        let simfunc: SimFunc = similarity.to_simfunc(true);
        let payloads_len = compute.core.payloads.len();
        let pb = crate::utils::simple_progressbar(payloads_len, "Computing similarity scores", multi);
        let mut payloads: Vec<SimScore> = Vec::with_capacity(payloads_len);

        #[rustfmt::skip]
        compute.core.payloads
            .into_par_iter()
            .progress_with(pb)
            .map(|sfvec| simfunc(&sfvec, &query))
            .collect_into_vec(&mut payloads);

        let mut roots = Vec::with_capacity(compute.core.roots.len());
        fn dive_recursive(
            store_node: Node<()>,
            above: &mut Vec<Node<BranchExtra>>,
            simscores: &[SimScore],
        ) -> Float {
            match store_node {
                Node::Branch(branch) => {
                    let mut current: Vec<Node<BranchExtra>> = Vec::with_capacity(branch.children.len());
                    let (n, sum) = branch
                        .children
                        .into_iter()
                        .map(|store_node| dive_recursive(store_node, &mut current, simscores))
                        .fold((0.0, 0.0), |(n, sum), x| (n + 1.0, sum + x));
                    let mean = sum / n;
                    above.push(Node::new_branch(
                        branch.name,
                        current.into_boxed_slice(),
                        BranchExtra { mean },
                    ));
                    mean
                }
                Node::Leaf(leaf) => {
                    let index = leaf.payload_id;
                    above.push(Node::Leaf(leaf));
                    unsafe { **simscores.get_unchecked(index) }
                }
            }
        }
        for store_node in compute.core.roots {
            dive_recursive(store_node, &mut roots, &payloads);
        }

        Ok(Self {
            core: TaxTreeCore {
                gene: compute.core.gene,
                highest_rank: compute.core.highest_rank,
                roots: roots.into_boxed_slice(),
                payloads: payloads.into_boxed_slice(),
            },
        })
    }

    pub fn cut(
        &self,
        nb_best: usize,
    ) -> Self {
        let best: SparseVec<SimScore> = self
            .core
            .gather_leaves()
            .iter()
            .map(|&leaf| (leaf.payload_id, unsafe { self.core.payloads.get_unchecked(leaf.payload_id) }))
            .sorted_by_key(|&(_, score)| score)
            .rev()
            .take(nb_best)
            .fold(SparseVec::new(usize::MAX), |mut acc, (pid, score)| {
                acc.insert(pid, *score);
                acc
            });

        fn dive_recursive(
            node: &Node<BranchExtra>,
            above: &mut Vec<Node<BranchExtra>>,
            new_payloads: &mut Vec<SimScore>,
            best: &SparseVec<SimScore>,
        ) {
            match node {
                Node::Branch(branch) => {
                    let mut add_self_flag: bool = false;
                    let mut current: Vec<Node<BranchExtra>> = Vec::new();
                    let mut leaves: Vec<&Leaf> = Vec::new();
                    for node in branch.children.iter() {
                        node.gather_leaves(&mut leaves);
                        if leaves.iter().any(|&leaf| best.keys().contains(&leaf.payload_id)) {
                            add_self_flag = true;
                            dive_recursive(node, &mut current, new_payloads, best);
                        }
                        leaves.clear();
                    }
                    if add_self_flag {
                        above.push(Node::new_branch(
                            branch.name.clone(),
                            current.into_boxed_slice(),
                            branch.extra,
                        ));
                    }
                }
                Node::Leaf(leaf) => {
                    above.push(Node::new_leaf(leaf.name.clone(), new_payloads.len()));
                    new_payloads.push(*best.get(leaf.payload_id).unwrap());
                }
            }
        }

        let mut new_payloads = Vec::new();
        let mut roots = Vec::new();
        for root in self.core.roots.iter() {
            dive_recursive(root, &mut roots, &mut new_payloads, &best);
        }

        TaxTreeResults {
            core: TaxTreeCore {
                gene: self.core.gene.clone(),
                highest_rank: self.core.highest_rank,
                roots: roots.into_boxed_slice(),
                payloads: new_payloads.into_boxed_slice(),
            },
        }
    }
}

impl std::fmt::Display for TaxTreeResults {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        writeln!(f, "{}", self.core)
    }
}
