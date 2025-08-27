#![allow(non_snake_case)]

// Imports
use std::{
    collections::HashMap,
    fs::OpenOptions,
    hint::unreachable_unchecked,
    io::Read,
    path::{Path, PathBuf},
    str::FromStr,
};

use bincode::{Decode, Encode};
use bio_seq::prelude::*;
#[cfg(feature = "indicatif")]
use indicatif::ParallelProgressIterator;
use itertools::Itertools;
use once_cell::unsync::OnceCell;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::{
    data::SFVec,
    tax::{Error, clade, clade::Rank, compute_sparse_kmer_counts_for_fasta_seq},
};

#[derive(Encode, Decode, PartialEq)]
pub enum StoreNode {
    Branch { name: Box<str>, children: Box<[StoreNode]> },
    Leaf { name: Box<str>, seq: Seq<Iupac>, pre_comp: Vec<SFVec> },
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
                write!(f, "({name}, seqlen={}, nb. pre-comp = {})", seq.len(), pre_comp.len())
            }
        }
    }
}

#[derive(Debug, Encode, Decode)]
pub struct TaxTreeStore {
    pub gene: String,
    pub filepath: PathBuf,
    pub highest_rank: Rank,
    pub pre_computed: Vec<usize>,
    pub roots: Box<[StoreNode]>,
}

impl TaxTreeStore {
    const BINCODE_CONFIG: bincode::config::Configuration = bincode::config::standard();
    pub const FILE_EXTENSION: &str = "myrtree";

    pub fn encode<W: std::io::Write>(
        &self,
        compression_level: i32,
        multithreading_flag: bool,
        nb_threads_available: usize,
        output: W,
    ) -> Result<W, Error> {
        #[cfg(feature = "indicatif")]
        let spinner = crate::utils::simple_spinner(
            Some(format!(
                "Encoding the '{}' gene taxonomic tree to '{}'",
                self.gene,
                self.filepath.display()
            )),
            Some(200),
        );

        let mut encoder = zstd::Encoder::new(output, compression_level)?;
        if multithreading_flag {
            encoder.multithread(nb_threads_available as u32)?;
        }
        bincode::encode_into_std_write(self, &mut encoder, Self::BINCODE_CONFIG)?;
        let output = encoder.finish()?;
        #[cfg(feature = "indicatif")]
        spinner.finish();
        Ok(output)
    }

    pub fn encode_to_file(
        &self,
        compression_level: i32,
        multithreading_flag: bool,
        nb_threads_available: usize,
    ) -> Result<(), Error> {
        let file = self.encode(
            compression_level,
            multithreading_flag,
            nb_threads_available,
            std::fs::OpenOptions::new().create(true).write(true).truncate(true).open(&self.filepath)?,
        )?;

        file.sync_all().map_err(Error::from)
    }

    pub fn decode<R: std::io::Read>(input: R) -> Result<Self, Error> {
        let mut decoder = zstd::Decoder::new(input)?;
        bincode::decode_from_std_read(&mut decoder, Self::BINCODE_CONFIG).map_err(Error::from)
    }

    pub fn decode_from_file<Q: AsRef<Path>>(filepath: Q) -> Result<Self, Error> {
        let filepath = filepath.as_ref().to_path_buf();
        let mut tree = Self::decode(std::fs::OpenOptions::new().read(true).open(&filepath)?)?;
        tree.filepath = filepath;
        Ok(tree)
    }

    /// Very heavy, do this only once if possible
    pub fn load_from_fasta_file<Q: AsRef<Path>>(
        input_filepath: Q,
        store_filepath: Option<Q>,
        gene: impl ToString,
        // the first element is the list of `k` for which to pre-compute sparse k-mer counts,
        // the second element is the `max_consecutive_N_before_cutoff` parameter
        pre_compute_kcounts: Option<(&[usize], usize)>,
    ) -> Result<Self, Error> {
        let store_filepath = if let Some(fp) = store_filepath {
            fp.as_ref().to_path_buf()
        } else {
            let mut fp = input_filepath.as_ref().to_path_buf();
            let file_name = fp.file_name().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();
            fp.set_file_name(file_name);
            fp.set_extension(Self::FILE_EXTENSION);
            fp
        };

        let mut file = OpenOptions::new().read(true).open(input_filepath)?;
        let mut fasta = match file.metadata() {
            Ok(metadata) => String::with_capacity(metadata.len() as usize),
            Err(e) => {
                eprintln!("Failed to get file metadata, {e}");
                String::new()
            }
        };
        file.read_to_string(&mut fasta)?;
        Self::load_from_fasta_string(fasta, store_filepath, gene, pre_compute_kcounts)
    }

    pub fn load_from_fasta_string(
        fasta: String,
        store_filepath: PathBuf,
        gene: impl ToString,
        // the first element is the list of `k` for which to pre-compute sparse k-mer counts,
        // the second element is the `max_consecutive_N_before_cutoff` parameter
        pre_compute_kcounts: Option<(&[usize], usize)>,
    ) -> Result<Self, Error> {
        /// Phylogenetic rank name stack (i.e., 'species' is at the top of the stack (end of vector), and 'domain' is at the bottom of the stack (start of vector))
        type Stack = Vec<Box<str>>;

        let mut lines = fasta.lines().enumerate().peekable();

        let highest_rank: OnceCell<Rank> = OnceCell::new();
        let mut leaves_and_stacks: Vec<(StoreNode, Stack)> = Vec::new();

        #[cfg(feature = "indicatif")]
        let spinner = crate::utils::simple_spinner(None, Some(200));

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
                Some((k_values, _)) => Vec::with_capacity(k_values.len()),
                None => Vec::with_capacity(0),
            };
            leaves_and_stacks.push((StoreNode::Leaf { name, seq, pre_comp }, stack));
            string_seq.clear();
        }

        #[cfg(feature = "indicatif")]
        spinner.finish();

        if let Some((k_values, max_consecutive_N_before_gap)) = pre_compute_kcounts {
            for &k in k_values.iter().unique() {
                #[cfg(feature = "indicatif")]
                let pbar_len = leaves_and_stacks.len();

                #[cfg(feature = "indicatif")]
                let parallel_iterator = leaves_and_stacks
                    .par_iter_mut()
                    .progress_with(crate::utils::simple_progressbar(pbar_len, format!("for k={k}")));

                #[cfg(not(feature = "indicatif"))]
                let parallel_iterator = leaves_and_stacks.par_iter_mut();

                parallel_iterator.for_each(|(leaf, _)| match leaf {
                    StoreNode::Leaf { name: _, seq, pre_comp } => {
                        pre_comp.push(compute_sparse_kmer_counts_for_fasta_seq(
                            seq,
                            k,
                            max_consecutive_N_before_gap,
                        ));
                    }
                    _ => unsafe { unreachable_unchecked() },
                });
            }
        }

        // Now we get to the fun part, start at the bottom of the ranks (at species) and move upwards while bundling together
        #[cfg(feature = "indicatif")]
        let spinner =
            crate::utils::simple_spinner(Some("Bundling everything together".to_string()), Some(200));

        let highest_rank =
            *highest_rank.get().expect("Totally invalid or empty file (highest rank was not set)");
        let mut nodes_and_stacks: Vec<(StoreNode, Stack)> = Vec::new();

        for curr_stack_len in (1..highest_rank as usize).rev() {
            let mut super_node_name_to_nodes_and_stack: HashMap<Box<str>, (Vec<StoreNode>, Stack)> =
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

        spinner.finish();

        Ok(Self {
            gene: gene.to_string(),
            filepath: store_filepath,
            highest_rank,
            roots,
            pre_computed: pre_compute_kcounts.map(|(k_values, _)| k_values.into()).unwrap_or_default(),
        })
    }

    /// Unsafe as the leaves must be `StoreNode::Leaf`
    pub(crate) unsafe fn compute_and_append_kmer_counts(
        mut leaves: Vec<&mut StoreNode>,
        k: usize,
        max_consecutive_N_before_gap: usize,
    ) {
        #[cfg(feature = "indicatif")]
        let pbar_len = leaves.len();

        #[cfg(feature = "indicatif")]
        let parallel_iterator = leaves
            .par_iter_mut()
            .progress_with(crate::utils::simple_progressbar(pbar_len, format!("for k={k}")));

        #[cfg(not(feature = "indicatif"))]
        let parallel_iterator = leaves_and_stacks.par_iter_mut();

        parallel_iterator.for_each(|leaf| match leaf {
            StoreNode::Leaf { name: _, seq, pre_comp } => {
                pre_comp.push(compute_sparse_kmer_counts_for_fasta_seq(seq, k, max_consecutive_N_before_gap));
            }
            _ => unsafe { unreachable_unchecked() },
        });
    }

    /// Unsafe as the leaves must be `StoreNode::Leaf`
    pub(crate) unsafe fn compute_and_overwrite_kmer_counts(
        mut leaves: Vec<&mut StoreNode>,
        k: usize,
        max_consecutive_N_before_gap: usize,
    ) {
        #[cfg(feature = "indicatif")]
        let pbar_len = leaves.len();

        #[cfg(feature = "indicatif")]
        let parallel_iterator = leaves
            .par_iter_mut()
            .progress_with(crate::utils::simple_progressbar(pbar_len, format!("for k={k}")));

        #[cfg(not(feature = "indicatif"))]
        let parallel_iterator = leaves_and_stacks.par_iter_mut();

        parallel_iterator.for_each(|leaf| match leaf {
            StoreNode::Leaf { name: _, seq, pre_comp } => {
                *pre_comp =
                    vec![compute_sparse_kmer_counts_for_fasta_seq(seq, k, max_consecutive_N_before_gap)]
            }
            _ => unsafe { unreachable_unchecked() },
        });
    }

    pub fn gather_leaves(&self) -> Vec<&StoreNode> {
        fn recursive_dive<'a>(
            node: &'a StoreNode,
            output: &mut Vec<&'a StoreNode>,
        ) {
            match node {
                StoreNode::Branch { name: _, children } => {
                    for node in children {
                        recursive_dive(node, output);
                    }
                }
                leaf @ StoreNode::Leaf { name: _, seq: _, pre_comp: _ } => output.push(leaf),
            }
        }
        let mut output = Vec::new();
        for root in self.roots.iter() {
            recursive_dive(root, &mut output);
        }
        output
    }

    pub fn gather_leaves_mut(&mut self) -> Vec<&mut StoreNode> {
        fn recursive_dive<'a>(
            node: &'a mut StoreNode,
            output: &mut Vec<&'a mut StoreNode>,
        ) {
            match node {
                StoreNode::Branch { name: _, children } => {
                    for node in children {
                        recursive_dive(node, output);
                    }
                }
                leaf @ StoreNode::Leaf { name: _, seq: _, pre_comp: _ } => output.push(leaf),
            }
        }
        let mut output = Vec::new();
        for root in self.roots.iter_mut() {
            recursive_dive(root, &mut output);
        }
        output
    }
}

#[cfg(test)]
mod test {
    use indoc::indoc;

    use super::*;

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

        let tree = TaxTreeStore::load_from_fasta_string(fasta, "".into(), "test", None).unwrap();

        let mut tree_leaves = tree.gather_leaves();
        tree_leaves.sort_by_key(|&node| match node {
            StoreNode::Leaf { name, seq: _, pre_comp: _ } => name,
            _ => unsafe { unreachable_unchecked() },
        });

        let expected_tree = TaxTreeStore {
            gene: "test".to_string(),
            filepath: PathBuf::new(),
            highest_rank: Rank::Order,
            pre_computed: vec![],
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
