// Imports
use crate::tax::clade::Rank;
use crate::tax::compute_sparse_kmer_counts_for_fasta_seq;
use crate::{data::SFVec, tax::clade};
use bincode::{Decode, Encode};
use bio_seq::prelude::*;
use itertools::Itertools;
use once_cell::unsync::OnceCell;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::hint::unreachable_unchecked;
use std::io::Read;
use std::{collections::HashMap, fs::OpenOptions, path::Path, str::FromStr};
use thiserror::Error;

#[cfg(feature = "indicatif")]
use indicatif::ParallelProgressIterator;

#[derive(Encode, Decode, PartialEq)]
pub enum StoreNode {
    Branch { name: Box<str>, children: Box<[StoreNode]> },
    Leaf { name: Box<str>, seq: Seq<Iupac>, pre_comp: Vec<(usize, SFVec)> },
}

impl core::fmt::Debug for StoreNode {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            Self::Branch { name, children } => {
                write!(f, "{name}: {children:?}")
            }
            Self::Leaf { name, seq, pre_comp } => {
                write!(f, "({name}, ...)")
            }
        }
    }
}

#[derive(Debug, Encode, Decode)]
pub struct TaxTreeStore {
    gene: String,
    roots: Box<[StoreNode]>,
    highest_rank: Rank,
}

impl TaxTreeStore {
    /// Very heavy, do this only once if possible
    pub fn load_from_fasta_file<Q: AsRef<Path>>(
        filepath: Q,
        gene: impl ToString,
        pre_compute_kcounts: Option<&[usize]>,
        pre_compute_kcounts_max_consecutive_N_before_cutoff: Option<usize>,
    ) -> Result<Self, Error> {
        let mut file = OpenOptions::new().read(true).open(filepath)?;
        let mut fasta = match file.metadata() {
            Ok(metadata) => String::with_capacity(metadata.len() as usize),
            Err(e) => {
                eprintln!("Failed to get file metadata, {e}");
                String::new()
            }
        };
        file.read_to_string(&mut fasta)?;
        Self::load_from_fasta_string(
            fasta,
            gene,
            pre_compute_kcounts,
            pre_compute_kcounts_max_consecutive_N_before_cutoff,
        )
    }

    pub fn load_from_fasta_string(
        fasta: String,
        gene: impl ToString,
        pre_compute_kcounts: Option<&[usize]>,
        pre_compute_kcounts_max_consecutive_N_before_cutoff: Option<usize>,
    ) -> Result<Self, Error> {
        let mut lines = fasta.lines().enumerate().peekable();

        let highest_rank: OnceCell<Rank> = OnceCell::new();
        let mut leaves_and_stacks: Vec<(StoreNode, Vec<Box<str>>)> = Vec::new();
        //let stack_size: OnceCell<usize> = OnceCell::new();

        #[cfg(feature = "indicatif")]
        let spinner = indicatif::ProgressBar::new_spinner().with_style(
            indicatif::ProgressStyle::default_spinner()
                .tick_chars("⊶⊷✔")
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );

        #[cfg(feature = "indicatif")]
        spinner.enable_steady_tick(std::time::Duration::from_millis(200));

        let mut string_seq = String::new();
        while let Some((lidx, text_line)) = lines.next() {
            #[cfg(feature = "indicatif")]
            spinner.set_message(format!("working on line n°{lidx}"));

            let (h_rank, _, mut stack) =
                clade::Parsed::from_str(text_line).map_err(|e| Error::CladeParse(e, lidx + 1))?.uncurl();

            let expected_highest_rank = *highest_rank.get_or_init(|| h_rank);
            if expected_highest_rank != h_rank {
                return Err(Error::HighestRankMismatch(expected_highest_rank, h_rank, lidx + 1));
            }

            while let Some((_, next_line)) = lines.peek() {
                if next_line.starts_with('>') {
                    break;
                }
                let dna_line = unsafe { lines.next().unwrap_unchecked().1 };
                string_seq.push_str(dna_line);
            }

            let name = unsafe { stack.pop().unwrap_unchecked() }; // Safe as an Ok(...) from `clade::parse_str` means the vector isn't empty
            let seq: Seq<Iupac> = Seq::from_str(&string_seq).map_err(|e| Error::BioSeq(e, lidx + 1))?;
            let pre_comp = match pre_compute_kcounts {
                Some(ks) => Vec::with_capacity(ks.len()),
                None => Vec::with_capacity(0),
            };
            leaves_and_stacks.push((StoreNode::Leaf { name, seq, pre_comp }, stack));
            string_seq.clear();
        }

        #[cfg(feature = "indicatif")]
        spinner.finish_with_message("Gathered tree leaves");

        if let Some(ks) = pre_compute_kcounts {
            #[allow(non_snake_case)]
            let max_consecutive_N_before_cutoff =
                pre_compute_kcounts_max_consecutive_N_before_cutoff.unwrap_or(4);
            for &k in ks.iter().unique() {
                #[cfg(feature = "indicatif")]
                let progress_bar = indicatif::ProgressBar::new(leaves_and_stacks.len() as u64)
                    .with_style(
                        indicatif::ProgressStyle::with_template(&format!(
                            "{{msg}} [{{elapsed_precise}}] {{bar:40.cyan/blue}} {{pos}} for k={k}"
                        ))
                        .unwrap(),
                    )
                    .with_message("⋆")
                    .with_finish(indicatif::ProgressFinish::WithMessage(crate::utils::greenify("✔").into()));

                #[cfg(feature = "indicatif")]
                let parallel_iterator = leaves_and_stacks.par_iter_mut().progress_with(progress_bar);

                #[cfg(not(feature = "indicatif"))]
                let parallel_iterator = leaves_and_stacks.par_iter_mut();

                parallel_iterator.for_each(|(leaf, _)| match leaf {
                    StoreNode::Leaf { name: _, seq, pre_comp } => {
                        pre_comp.push((
                            k,
                            compute_sparse_kmer_counts_for_fasta_seq(seq, k, max_consecutive_N_before_cutoff),
                        ));
                    }
                    _ => unsafe { unreachable_unchecked() },
                });
            }
        }

        // Now we get to the fun part, start at the bottom of the ranks (at species) and move upwards while bundling together
        let highest_rank =
            *highest_rank.get().expect("Totally invalid or empty file (highest rank was not set)");

        let mut nodes_and_stacks: Vec<(StoreNode, Vec<Box<str>>)> = Vec::new();

        for curr_stack_len in (1..highest_rank as usize).rev() {
            let mut super_node_name_to_nodes_and_stack: HashMap<Box<str>, (Vec<StoreNode>, Vec<Box<str>>)> =
                HashMap::new();

            let mut remaining_leaves_and_stacks = Vec::new();
            leaves_and_stacks
                .into_iter()
                .filter_map(|las| {
                    if las.1.len() == curr_stack_len {
                        Some(las)
                    } else {
                        remaining_leaves_and_stacks.push(las);
                        None
                    }
                })
                .chain(nodes_and_stacks.into_iter())
                .for_each(|(node, mut stack)| {
                    let super_name = unsafe { stack.pop().unwrap_unchecked() };
                    super_node_name_to_nodes_and_stack
                        .entry(super_name)
                        .or_insert((Vec::new(), stack))
                        .0
                        .push(node);
                });

            leaves_and_stacks = remaining_leaves_and_stacks;

            nodes_and_stacks = super_node_name_to_nodes_and_stack
                .into_iter()
                .map(|(name, (nodes, stack))| {
                    (StoreNode::Branch { name, children: nodes.into_boxed_slice() }, stack)
                })
                .collect_vec();
        }

        let roots = nodes_and_stacks
            .into_iter()
            .chain(leaves_and_stacks)
            .map(|(node, _)| node)
            .collect_vec()
            .into_boxed_slice();

        //println!("{:?}", roots.first().unwrap());

        Ok(Self { gene: gene.to_string(), roots, highest_rank })
    }

    pub fn gather_leaves(&self) -> Vec<&StoreNode> {
        fn recursive_gather_leaves<'a>(
            node: &'a StoreNode,
            output: &mut Vec<&'a StoreNode>,
        ) {
            match node {
                StoreNode::Branch { name: _, children } => {
                    for node in children {
                        recursive_gather_leaves(node, output);
                    }
                }
                leaf @ StoreNode::Leaf { name: _, seq: _, pre_comp: _ } => output.push(leaf),
            }
        }
        let mut output = Vec::new();
        for root in self.roots.iter() {
            recursive_gather_leaves(root, &mut output);
        }
        output
    }
}

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    IO(#[from] std::io::Error),
    #[error("Failed to parse taxonomic identity of the record starting on line {1}, {0}")]
    CladeParse(clade::ParsingError, usize),
    #[error("Bio-seq parsing error for the record starting on line {1}, {0}")]
    BioSeq(ParseBioError, usize),
    #[error("Expected the highest rank to be {0}, got {1} instead for the record starting on line {2}")]
    HighestRankMismatch(Rank, Rank, usize),
}

#[cfg(test)]
mod test {
    use super::*;
    use indoc::indoc;

    #[test]
    fn load_fasta_test() {
        let fasta = indoc! {"
            >tax={ o:order_1; f:fam_1, g:genus_1, s:mouse }
            ACTG
            ACTG
            >tax={ o:order_1; f:fam_1, g:genus_1, s:rat }
            A
            >tax={ o:order_1; f:fam_1, g:genus_1 }
            A
            >tax={ o: order_2 }
            A
            >tax={ o: order_3, f:fam_2, g:genus_2, s: ant }
            A
            >tax={ o: order_3, f:fam_2, g:genus_2, s: termite }
            A
        "}
        .to_string();

        let tree = TaxTreeStore::load_from_fasta_string(fasta, "test", None, None).unwrap();

        let mut tree_leaves = tree.gather_leaves();
        tree_leaves.sort_by_key(|&node| match node {
            StoreNode::Leaf { name, seq: _, pre_comp: _ } => name,
            _ => unsafe { unreachable_unchecked() },
        });

        let expected_tree = TaxTreeStore {
            gene: "test".to_string(),
            roots: Box::from([
                StoreNode::Branch {
                    name: Box::from("order_1"),
                    children: Box::from([StoreNode::Branch {
                        name: Box::from("fam_1"),
                        children: Box::from([
                            StoreNode::Branch {
                                name: Box::from("genus_1"),
                                children: Box::from([
                                    StoreNode::Leaf {
                                        name: Box::from("mouse"),
                                        seq: iupac!("ACTGACTG").to_owned(),
                                        pre_comp: vec![],
                                    },
                                    StoreNode::Leaf {
                                        name: Box::from("rat"),
                                        seq: iupac!("A").to_owned(),
                                        pre_comp: vec![],
                                    },
                                ]),
                            },
                            StoreNode::Leaf {
                                name: Box::from("genus_1"),
                                seq: iupac!("A").to_owned(),
                                pre_comp: vec![],
                            },
                        ]),
                    }]),
                },
                StoreNode::Leaf { name: Box::from("order_2"), seq: iupac!("A").to_owned(), pre_comp: vec![] },
                StoreNode::Branch {
                    name: Box::from("order_3"),
                    children: Box::from([StoreNode::Branch {
                        name: Box::from("fam_2"),
                        children: Box::from([StoreNode::Branch {
                            name: Box::from("genus_2"),
                            children: Box::from([
                                StoreNode::Leaf {
                                    name: Box::from("ant"),
                                    seq: iupac!("A").to_owned(),
                                    pre_comp: vec![],
                                },
                                StoreNode::Leaf {
                                    name: Box::from("termite"),
                                    seq: iupac!("A").to_owned(),
                                    pre_comp: vec![],
                                },
                            ]),
                        }]),
                    }]),
                },
            ]),
            highest_rank: Rank::Order,
        };

        let mut expected_tree_leaves = expected_tree.gather_leaves();
        expected_tree_leaves.sort_by_key(|&node| match node {
            StoreNode::Leaf { name, seq: _, pre_comp: _ } => name,
            _ => unsafe { unreachable_unchecked() },
        });

        println!("{tree_leaves:?}\n{expected_tree_leaves:?}");

        assert_eq!(tree_leaves, expected_tree_leaves);
    }
}
