#![allow(non_snake_case)]

// Imports
use std::{
    collections::HashMap,
    fs::OpenOptions,
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
    tax::{
        Error,
        clade::{self, Rank},
        compute_sparse_kmer_counts_for_fasta_seq,
        core::{Leaf, Node, TaxTreeCore},
    },
};

#[derive(Encode, Decode)]
pub struct StorePayload {
    pub(crate) seq: Seq<Iupac>,
    pub(crate) pre_comp: Vec<SFVec>,
}

impl StorePayload {
    fn new(
        seq: Seq<Iupac>,
        pre_comp: Vec<SFVec>,
    ) -> Self {
        Self { seq, pre_comp }
    }
}

#[derive(Encode, Decode)]
pub struct TaxTreeStore {
    pub(crate) core: TaxTreeCore<(), StorePayload>,
    pub(crate) filepath: PathBuf,
    pub(crate) pre_computed: Vec<usize>,
}

impl TaxTreeStore {
    const BINCODE_CONFIG: bincode::config::Configuration = bincode::config::standard();
    pub const FILE_EXTENSION: &str = "myrtree";

    pub fn encode<W: std::io::Write>(
        &self,
        zstd_compression_level: i32,
        zstd_multithreading_opt: Option<u32>,
        output: W,
    ) -> Result<W, Error> {
        let mut encoder = zstd::Encoder::new(output, zstd_compression_level)?;
        if let Some(n_workers) = zstd_multithreading_opt {
            encoder.multithread(n_workers)?;
        }
        bincode::encode_into_std_write(self, &mut encoder, Self::BINCODE_CONFIG)?;
        encoder.finish().map_err(Error::from)
    }

    pub fn encode_to_file(
        &self,
        zstd_compression_level: i32,
        zstd_multithreading_opt: Option<u32>,
    ) -> Result<(), Error> {
        #[cfg(feature = "indicatif")]
        let spinner = crate::utils::simple_spinner(
            Some(format!(
                "Encoding the '{}' gene taxonomic tree to '{}'",
                self.core.gene,
                self.filepath.display()
            )),
            Some(200),
        );
        let file = self.encode(
            zstd_compression_level,
            zstd_multithreading_opt,
            std::fs::OpenOptions::new().create(true).write(true).truncate(true).open(&self.filepath)?,
        )?;
        file.sync_all()?;

        #[cfg(feature = "indicatif")]
        spinner.finish();

        Ok(())
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
        type StoreNode = Node<()>;

        let mut lines = fasta.lines().enumerate().peekable();

        let highest_rank: OnceCell<Rank> = OnceCell::new();
        let mut leaves_and_stacks: Vec<(StoreNode, Stack)> = Vec::new();
        let mut payloads: Vec<StorePayload> = Vec::new();
        let mut payload_id: usize = 0;

        let pre_comp = match pre_compute_kcounts {
            Some((k_values, _)) => Vec::with_capacity(k_values.len()),
            None => Vec::with_capacity(0),
        };

        #[cfg(feature = "indicatif")]
        let spinner = crate::utils::simple_spinner(None, Some(200));

        let mut string_seq = String::new();

        'outer: while let Some((lidx, text_line)) = lines.next() {
            #[cfg(feature = "indicatif")]
            spinner.set_message(format!("working on line n°{lidx}"));

            let (h_rank, _, mut stack) =
                match clade::Parsed::from_str(text_line).map_err(|e| Error::CladeParse(e, lidx + 1)) {
                    Ok(parsed) => parsed.uncurl(),
                    Err(e) => {
                        eprintln!("{e}");
                        while let Some((_, next_line)) = lines.peek() {
                            if next_line.starts_with('>') {
                                continue 'outer;
                            }
                            lines.next();
                        }
                        break;
                    }
                };

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

            leaves_and_stacks.push((StoreNode::new_leaf(name, payload_id), stack));
            payloads.push(StorePayload::new(seq, pre_comp.clone()));
            payload_id += 1;
            string_seq.clear();
        }

        #[cfg(feature = "indicatif")]
        spinner.finish();

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
                    (StoreNode::new_branch(name, nodes.into_boxed_slice(), ()), stack)
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

        if let Some((k_values, max_consecutive_N_before_gap)) = pre_compute_kcounts {
            let mut tax_tree_store = Self {
                core: TaxTreeCore::new(gene.to_string(), highest_rank, roots, payloads.into_boxed_slice()),
                filepath: store_filepath,
                pre_computed: pre_compute_kcounts.map(|(k_values, _)| k_values.into()).unwrap_or_default(),
            };

            for &k in k_values {
                tax_tree_store.compute_and_append_kmer_counts(k, max_consecutive_N_before_gap);
            }

            Ok(tax_tree_store)
        } else {
            Ok(Self {
                core: TaxTreeCore::new(gene.to_string(), highest_rank, roots, payloads.into_boxed_slice()),
                filepath: store_filepath,
                pre_computed: pre_compute_kcounts.map(|(k_values, _)| k_values.into()).unwrap_or_default(),
            })
        }
    }

    pub(crate) fn compute_and_append_kmer_counts(
        &mut self,
        k: usize,
        max_consecutive_N_before_gap: usize,
    ) {
        #[cfg(feature = "indicatif")]
        let pbar_len = self.core.payloads.len();

        #[cfg(feature = "indicatif")]
        let parallel_iterator = self
            .core
            .payloads
            .par_iter_mut()
            .progress_with(crate::utils::simple_progressbar(pbar_len, format!("for k={k}")));

        #[cfg(not(feature = "indicatif"))]
        let parallel_iterator = self.core.payloads.par_iter_mut();

        parallel_iterator.for_each(|sp| {
            let mut kmer_freqs =
                compute_sparse_kmer_counts_for_fasta_seq(&sp.seq, k, max_consecutive_N_before_gap);
            kmer_freqs /= kmer_freqs.sum();
            sp.pre_comp.push(kmer_freqs);
        });
    }

    pub(crate) fn compute_and_overwrite_kmer_counts(
        &mut self,
        k: usize,
        max_consecutive_N_before_gap: usize,
    ) {
        #[cfg(feature = "indicatif")]
        let pbar_len = self.core.payloads.len();

        #[cfg(feature = "indicatif")]
        let parallel_iterator = self
            .core
            .payloads
            .par_iter_mut()
            .progress_with(crate::utils::simple_progressbar(pbar_len, format!("for k={k}")));

        #[cfg(not(feature = "indicatif"))]
        let parallel_iterator = self.core.payloads.par_iter_mut();

        parallel_iterator.for_each(|sp| {
            let mut kmer_freqs =
                compute_sparse_kmer_counts_for_fasta_seq(&sp.seq, k, max_consecutive_N_before_gap);
            kmer_freqs /= kmer_freqs.sum();
            let pre_comp = &mut sp.pre_comp;
            *pre_comp = vec![kmer_freqs]
        });
    }
}
