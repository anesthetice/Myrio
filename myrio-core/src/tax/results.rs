#![allow(non_snake_case)]

use indicatif::MultiProgress;
use itertools::Itertools;
// Imports
use indicatif::ParallelProgressIterator;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelIterator;
use rayon::iter::ParallelIterator;

use crate::similarity::SimFunc;
use crate::tax::compute::TaxTreeCompute;
use crate::tax::core::Node;
use crate::{
    data::SFVec,
    similarity::SimScore,
    tax::{Error, core::TaxTreeCore},
};

pub struct BranchExtra {
    mean: f64,
}

pub struct TaxTreeResults {
    core: TaxTreeCore<BranchExtra, SimScore>,
}

impl TaxTreeResults {
    pub fn from_compute_tree(
        query: SFVec,
        compute: TaxTreeCompute,
        simfunc: SimFunc,
        multi: Option<&MultiProgress>,
    ) -> Result<Self, Error> {
        let payloads_len = compute.core.payloads.len();
        let pb = crate::utils::simple_progressbar(payloads_len, "Computing similarity scores", multi);
        let mut payloads: Vec<SimScore> = Vec::with_capacity(payloads_len);
        #[rustfmt::skip]
        compute.core.payloads
            .into_par_iter()
            .progress_with(pb)
            .map(|sfvec| simfunc(&query, &sfvec))
            .collect_into_vec(&mut payloads);

        let mut roots = Vec::with_capacity(compute.core.roots.len());
        fn dive_recursive(
            store_node: Node<()>,
            above: &mut Vec<Node<BranchExtra>>,
            simscores: &[SimScore],
        ) -> f64 {
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

    pub fn test(&self) {
        let best = self
            .core
            .gather_leaves()
            .iter()
            .map(|leaf| (&leaf.name, self.core.payloads[leaf.payload_id]))
            .max_by_key(|(_, score)| *score)
            .unwrap();

        println!("for {}, best = {}, with a score of {}", self.core.gene, best.0, best.1);
    }
}
