// Imports
use std::ops::Neg;

use indicatif::{MultiProgress, ParallelProgressIterator};
use indoc::printdoc;
use itertools::Itertools;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator,
    ParallelIterator,
};

use crate::{
    data::{Float, SFVec, SparseVec},
    similarity::{SimFunc, SimScore, Similarity},
    tax::{
        Error,
        clade::Rank,
        compute::TaxTreeCompute,
        core::{Leaf, Node, TaxTreeCore},
    },
};

#[derive(Clone, Copy)]
pub struct BRes {
    pub nb_leaves: usize,
    pub mean: Option<Float>,
    pub max: Option<Float>,
    pub pool_score: Option<Float>,
}

impl std::fmt::Display for BRes {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        let mut s: String = String::new();
        s += &format!("🮰: {}", self.nb_leaves);
        if let Some(mean) = self.mean {
            s += &format!(", μ={mean:.3}")
        }
        if let Some(max) = self.max {
            s += &format!(", max={max:.3}")
        }
        if let Some(pool_score) = self.pool_score {
            s += &format!(", s={pool_score:.3}")
        }
        f.write_str(s.as_str())
    }
}

pub struct TaxTreeResults {
    pub core: TaxTreeCore<BRes, SimScore>,
}

impl TaxTreeResults {
    pub fn from_compute_tree(
        query: SFVec,
        ttcompute: TaxTreeCompute,
        similarity: Similarity,
        multi: Option<&MultiProgress>,
    ) -> Result<Self, Error> {
        let simfunc: SimFunc = similarity.to_simfunc(true);

        let payloads_len = ttcompute.core.payloads.len();
        let mut payloads: Vec<SimScore> = Vec::with_capacity(payloads_len);
        let pb = crate::utils::simple_progressbar(payloads_len, "Computing similarity scores", multi);

        #[rustfmt::skip]
        ttcompute.core.payloads
            .into_par_iter()
            .progress_with(pb)
            .map(|sfvec| simfunc(&sfvec, &query))
            .collect_into_vec(&mut payloads);

        let mut roots = Vec::with_capacity(ttcompute.core.roots.len());

        // This recursive function serves to convert the roots of TaxTreeCore<(), SFVec> used by TaxTreeCompute into the roots of TaxTreeCore<BRes, SimScore> used by TaxTreeResults
        fn dive_recursive(
            store_node: Node<()>,
            above: &mut Vec<Node<BRes>>,
        ) -> usize {
            match store_node {
                Node::Branch(store_branch) => {
                    let mut results_branch_children: Vec<Node<BRes>> =
                        Vec::with_capacity(store_branch.children.len());
                    let nb_leaves = store_branch
                        .children
                        .into_iter()
                        .map(|store_node| dive_recursive(store_node, &mut results_branch_children))
                        .sum::<usize>();
                    let results_branch = Node::new_branch(
                        store_branch.name,
                        results_branch_children.into_boxed_slice(),
                        BRes { nb_leaves, mean: None, max: None, pool_score: None },
                    );
                    above.push(results_branch);
                    nb_leaves
                }
                Node::Leaf(leaf) => {
                    above.push(Node::Leaf(leaf));
                    1_usize
                }
            }
        }
        for store_node in ttcompute.core.roots {
            dive_recursive(store_node, &mut roots);
        }

        // The incomplete results tree, notice that the payloads haven't been inserted yet (for borrow checker related reasons)
        let mut ttres = Self {
            core: TaxTreeCore {
                gene: ttcompute.core.gene,
                highest_rank: ttcompute.core.highest_rank,
                roots: roots.into_boxed_slice(),
                payloads: Box::from([]),
            },
        };

        let mut genus_vec = ttres.core.gather_branches_mut_at_rank(Rank::Genus);
        genus_vec.par_iter_mut().for_each(|genus| {
            let scores = genus
                .gather_leaves()
                .into_iter()
                .map(|leaf| unsafe { **payloads.get_unchecked(leaf.payload_id) })
                .collect_vec();

            genus.extra.mean.replace(scores.iter().sum::<Float>() / scores.len() as Float);
            genus.extra.max.replace(scores.iter().copied().reduce(Float::max).unwrap());

            // Softmax Pooling multiplied by a low `n` penalty
            static LAMBDA: Float = 1.5;
            static BETA: Float = 1.2;
            let pool_score = {
                let numerator = scores.iter().copied().map(|s| s * (LAMBDA * s).exp()).sum::<Float>();
                let denominator = scores.iter().copied().map(|s| (LAMBDA * s).exp()).sum::<Float>();
                let small_size_penalty_multiplier = scores.len() as Float / (scores.len() as Float + BETA);

                (numerator / denominator) * small_size_penalty_multiplier
            };
            genus.extra.pool_score.replace(pool_score);
        });

        {
            let n = genus_vec.len() as Float;
            let pool_scores = genus_vec
                .iter()
                .map(|branch| unsafe { branch.extra.pool_score.unwrap_unchecked() })
                .collect_vec();

            let pool_score_mean = pool_scores.iter().sum::<Float>() / n;
            #[rustfmt::skip]
            let pool_score_std =
                (n.powi(-1) * pool_scores.into_iter().map(|x| (x - pool_score_mean).powi(2)).sum::<Float>()).sqrt();

            static GAMMA: Float = 2.5;
            static DELTA: Float = 0.9;

            let first_genus = {
                let pos = genus_vec
                    .iter()
                    .position_max_by_key(|branch| unsafe {
                        SimScore::new_unchecked(branch.extra.pool_score.unwrap_unchecked())
                    })
                    .unwrap();
                genus_vec.remove(pos)
            };

            let second_genus = {
                let pos = genus_vec
                    .iter()
                    .position_max_by_key(|branch| unsafe {
                        SimScore::new_unchecked(branch.extra.pool_score.unwrap_unchecked())
                    })
                    .unwrap();
                genus_vec.remove(pos)
            };

            let x = (first_genus.extra.pool_score.unwrap() - second_genus.extra.pool_score.unwrap())
                / pool_score_std;

            let conf =
                first_genus.extra.pool_score.unwrap().powf(0.7) * x.powf(GAMMA) / (DELTA + x.powf(GAMMA));

            printdoc! {"
                1st genus: {} -> {:.3}
                2nd genus: {} -> {:.3}
                mean: {pool_score_mean:.3}, std: {pool_score_std:.3}, x={x:.3}
                => conf = {conf:.3}
            ",
            first_genus.name, first_genus.extra.pool_score.unwrap(),
            second_genus.name, second_genus.extra.pool_score.unwrap()
            };
        }

        // insert payloads before returning
        ttres.core.payloads = payloads.into_boxed_slice();

        Ok(ttres)
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
            .sorted_by_key(|&(_, score)| unsafe { SimScore::new_unchecked(score.neg()) })
            .take(nb_best)
            .fold(SparseVec::new(usize::MAX), |mut acc, (pid, score)| {
                acc.insert(pid, *score);
                acc
            });

        fn dive_recursive(
            node: &Node<BRes>,
            above: &mut Vec<Node<BRes>>,
            new_payloads: &mut Vec<SimScore>,
            best: &SparseVec<SimScore>,
        ) {
            match node {
                Node::Branch(branch) => {
                    let mut add_self_flag: bool = false;
                    let mut current: Vec<Node<BRes>> = Vec::new();
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
