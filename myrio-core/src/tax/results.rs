// Imports
use std::ops::Neg;

use console::{StyledObject, style};
use indicatif::{MultiProgress, ParallelProgressIterator};
use indoc::printdoc;
use itertools::Itertools;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator,
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

#[derive(Clone)]
pub struct BRes {
    pub nb_leaves: usize,
    pub mean: Float,
    pub max: Float,
    pub pool_score: Float,
    pub indicator: Option<String>,
}

impl std::fmt::Display for BRes {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        if let Some(indicator) = &self.indicator {
            write!(
                f,
                "{indicator} | ({}, μ={:.3}, max={:.3}, s={:.3})",
                self.nb_leaves, self.mean, self.max, self.pool_score
            )
        } else {
            write!(
                f,
                "({}, μ={:.3}, max={:.3}, s={:.3})",
                self.nb_leaves, self.mean, self.max, self.pool_score
            )
        }
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
        lambda_leaf: Float,
        lambda_branch: Float,
        mu: Float,
        gamma: Float,
        delta: Float,
        epsilon: Float,
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

        // This recursive function serves to convert the roots of TaxTreeCore<(), SFVec> used by TaxTreeCompute into the roots of TaxTreeCore<BRes, SimScore> used by TaxTreeResults. The BRes parameters `nb_leaves`, `mean`, `max` are all computed here, `pool_score` is set to the default of 0.0 and will be computed later
        fn dive_recursive(
            store_node: Node<()>,
            above: &mut Vec<Node<BRes>>,
            simscores: &[SimScore],
        ) -> (usize, Float, Float) // (number of leaves, sum, max)
        {
            match store_node {
                Node::Branch(store_branch) => {
                    let mut results_branch_children: Vec<Node<BRes>> =
                        Vec::with_capacity(store_branch.children.len());

                    let mut nb_leaves: usize = 0;
                    let mut sum: Float = 0.0;
                    let mut max: Float = 0.0;

                    for store_node in store_branch.children {
                        let (nb_leaves_i, sum_i, max_i) =
                            dive_recursive(store_node, &mut results_branch_children, simscores);

                        nb_leaves += nb_leaves_i;
                        sum += sum_i;
                        max = max.max(max_i);
                    }

                    let results_branch = Node::new_branch(
                        store_branch.name,
                        results_branch_children.into_boxed_slice(),
                        BRes {
                            nb_leaves,
                            mean: sum / nb_leaves as Float,
                            max,
                            pool_score: 0.0,
                            indicator: None,
                        },
                    );
                    above.push(results_branch);
                    (nb_leaves, sum, max)
                }
                Node::Leaf(leaf) => {
                    let score = **unsafe { simscores.get_unchecked(leaf.payload_id) };
                    above.push(Node::Leaf(leaf));
                    (1_usize, score, score)
                }
            }
        }
        for store_node in ttcompute.core.roots {
            dive_recursive(store_node, &mut roots, &payloads);
        }

        let ttres = Self {
            core: TaxTreeCore {
                gene: ttcompute.core.gene,
                highest_rank: ttcompute.core.highest_rank,
                roots: roots.into_boxed_slice(),
                payloads: payloads.into_boxed_slice(),
            },
        };
        let mut ttres = ttres.cut(100);

        // We extract the payloads to satisfy the borrow-checker, these will be inserted later at the end
        let payloads = std::mem::replace(&mut ttres.core.payloads, Box::from([]));

        fn softmax_pooling(
            scores: Vec<Float>,
            lambda: Float,
        ) -> Float {
            let mut numerator: Float = 0.0;
            let mut denominator: Float = 0.0;

            for s in scores.into_iter() {
                numerator += s * (lambda * s).exp();
                denominator += (lambda * s).exp();
            }

            numerator / denominator
        }

        for rank in Rank::collect_range_inclusive(Rank::Genus, ttres.core.highest_rank) {
            let mut rank_branch_vec = ttres.core.gather_branches_mut_at_rank(rank);

            rank_branch_vec.par_iter_mut().for_each(|rank_branch| {
                let mut sim_scores_below: Vec<Float> = Vec::new();
                let mut pool_scores_below: Vec<Float> = Vec::new();

                for sub_node in rank_branch.children.iter() {
                    match sub_node {
                        Node::Branch(sub_branch) => {
                            pool_scores_below.push(sub_branch.extra.pool_score);
                        }
                        Node::Leaf(leaf) => {
                            sim_scores_below.push(unsafe { **payloads.get_unchecked(leaf.payload_id) })
                        }
                    }
                }
                if !sim_scores_below.is_empty() {
                    // Softmax Pooling multiplied by a low `n` penalty
                    let n: Float = sim_scores_below.len() as Float;
                    let penalty = n / (n + mu);
                    pool_scores_below.push(softmax_pooling(sim_scores_below, lambda_leaf) * penalty);
                }
                rank_branch.extra.pool_score = softmax_pooling(pool_scores_below, lambda_branch);
            });

            if rank_branch_vec.is_empty() {
                continue;
            } else if rank_branch_vec.len() == 1 {
                let single_branch = rank_branch_vec.pop().unwrap();
                single_branch.extra.indicator = Some(style("◉: ∅").green().to_string());
                continue;
            }

            let n = rank_branch_vec.len() as Float;
            let pool_scores = rank_branch_vec.iter().map(|branch| branch.extra.pool_score).collect_vec();

            let pool_score_mean = pool_scores.iter().sum::<Float>() / n;
            #[rustfmt::skip]
            let pool_score_std =
                (n.powi(-1) * pool_scores.into_iter().map(|x| (x - pool_score_mean).powi(2)).sum::<Float>()).sqrt();

            let first = rank_branch_vec.remove(
                rank_branch_vec
                    .iter()
                    .position_max_by_key(|branch| SimScore::try_new(branch.extra.pool_score).unwrap())
                    .unwrap(),
            );

            let second = rank_branch_vec.remove(
                rank_branch_vec
                    .iter()
                    .position_max_by_key(|branch| SimScore::try_new(branch.extra.pool_score).unwrap())
                    .unwrap(),
            );

            let x = (first.extra.pool_score - second.extra.pool_score) / pool_score_std;

            let conf = first.extra.pool_score.powf(epsilon) * x.powf(gamma) / (delta + x.powf(gamma));

            let colorizer = match conf {
                0.0..0.25 => StyledObject::<String>::red,
                0.25..0.5 => StyledObject::<String>::yellow,
                0.5..0.75 => StyledObject::<String>::green,
                0.75..1.0 => StyledObject::<String>::cyan,
                _ => {
                    eprintln!("Warning: got unexpected confidence score");
                    StyledObject::magenta
                }
            };
            first.extra.indicator = Some(colorizer(style(format!("◉: {conf:.3}"))).to_string());

            #[cfg(debug_assertions)]
            printdoc! {"
                ## {rank}
                1st : {} -> {:.3} | 2nd : {} -> {:.3}
                mean: {pool_score_mean:.3}, std: {pool_score_std:.3}, x={x:.3}
                => conf = {conf:.3}\n
            ",
            first.name, first.extra.pool_score,
            second.name, second.extra.pool_score,
            };
        }

        // insert payloads before returning
        ttres.core.payloads = payloads;

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
                            branch.extra.clone(),
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
