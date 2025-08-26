// Imports
use crate::{
    data::SFVec,
    tax::{
        clade::Rank,
        store::{StoreNode, TaxTreeStore},
    },
};

#[derive(Debug)]
pub struct Branch {
    pub name: Box<str>,
    pub children: Box<[Node]>,
}

#[derive(Debug)]
pub struct Leaf {
    pub name: Box<str>,
    pub kmer_counts: SFVec,
}

#[derive(Debug)]
pub enum Node {
    Branch(Branch),
    Leaf(Leaf),
}

impl Node {
    pub fn new_branch(
        name: Box<str>,
        children: Box<[Node]>,
    ) -> Self {
        Self::Branch(Branch { name, children })
    }

    pub fn new_leaf(
        name: Box<str>,
        kmer_counts: SFVec,
    ) -> Self {
        Self::Leaf(Leaf { name, kmer_counts })
    }
}

#[derive(Debug)]
pub struct TaxTree {
    gene: String,
    k: usize,
    highest_rank: Rank,
    roots: Box<[Node]>,
}

impl TaxTree {
    pub fn from_store(
        mut store: TaxTreeStore,
        k: usize,
        cache_results: bool,
    ) -> Self {
        let sfvec_idx = store.pre_computed.iter().copied().find(|&k_| k == k_).unwrap_or_else(|| {
            let store_leaves_mut = store.gather_leaves_mut();
            if cache_results {
                unsafe { TaxTreeStore::compute_and_append_kmer_counts(store_leaves_mut, k, None) };
                store.encode_to_file(17);
                store.pre_computed.len()
            } else {
                unsafe { TaxTreeStore::compute_and_overwrite_kmer_counts(store_leaves_mut, k, None) };
                0
            }
        });

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

        let mut roots: Vec<Node> = Vec::with_capacity(store.roots.len());
        for store_root in store.roots.into_iter() {
            dive_recursive(store_root, &mut roots, sfvec_idx);
        }

        Self { gene: store.gene, k, highest_rank: store.highest_rank, roots: roots.into_boxed_slice() }
    }

    pub fn gather_leaves(&self) -> Vec<&Leaf> {
        fn recursive_dive<'a>(
            node: &'a Node,
            output: &mut Vec<&'a Leaf>,
        ) {
            match node {
                Node::Branch(branch) => {
                    for node in &branch.children {
                        recursive_dive(node, output);
                    }
                }
                Node::Leaf(leaf) => output.push(leaf),
            }
        }
        let mut output = Vec::new();
        for root in self.roots.iter() {
            recursive_dive(root, &mut output);
        }
        output
    }
}
