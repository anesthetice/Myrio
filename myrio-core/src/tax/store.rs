// Imports
use std::{
    collections::HashMap,
    fs::OpenOptions,
    io::{Read, Write},
    path::{Path, PathBuf},
    str::FromStr,
};

use bincode::{Decode, Encode};
use bio_seq::prelude::*;
use indicatif::{MultiProgress, ParallelProgressIterator};
use itertools::Itertools;
use once_cell::unsync::OnceCell;
use rand::{SeedableRng, rngs::SmallRng};
use rayon::{
    iter::{IndexedParallelIterator, IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator},
    slice::ParallelSlice,
};

use crate::{
    data::{Float, SparseVec},
    tax::{
        Error,
        clade::{self, Rank},
        compute_kmer_store_counts_for_fasta_seq,
        core::{Node, TaxTreeCore},
    },
};

#[derive(Debug, Clone, PartialEq, Encode, Decode)]
pub struct StorePayload {
    pub(crate) seq: Seq<Iupac>,
    pub(crate) kmer_store_counts_vec: Vec<(SparseVec<u16>, Float)>,
}

impl StorePayload {
    fn new(
        seq: Seq<Iupac>,
        kmer_store_counts_vec: Vec<(SparseVec<u16>, Float)>,
    ) -> Self {
        Self { seq, kmer_store_counts_vec }
    }
}

#[derive(Encode, Decode)]
pub struct TaxTreeStore {
    pub(crate) core: TaxTreeCore<(), StorePayload>,
    pub(crate) filepath: PathBuf,
    pub k_precomputed: Vec<usize>,
}

#[derive(Debug, Clone, Encode, Decode)]
pub struct TTSHeader {
    core_size: usize,
    payloads_chunk_sizes: Vec<usize>,
    k_precomputed: Vec<usize>,
}

impl TaxTreeStore {
    const BINCODE_CONFIG: bincode::config::Configuration = bincode::config::standard();
    pub const FILE_EXTENSION: &str = "myrtree";
    // "MYRIO-Ψ"
    const MAGIC_NUMBER: [u8; 8] = [77, 89, 82, 73, 79, 45, 206, 168];

    /// Self needs to be mutable so we can temporarily swap out the payloads
    pub fn to_bytes(
        &mut self,
        zstd_compression_level: i32,
    ) -> Result<Vec<u8>, Error> {
        // First we extract the payloads by swapping them with an empty boxed array
        let extracted_payloads = std::mem::replace(&mut self.core.payloads, Box::new([]));

        // Then we encode and compress the core (without any payloads)
        let compressed_encoded_core = zstd::encode_all(
            bincode::encode_to_vec(&self.core, Self::BINCODE_CONFIG)?.as_slice(),
            zstd_compression_level,
        )?;

        // Now we encode the payloads into chunks in parallel
        // We want roughly eight or so chunks, might make this editable in the future
        let chunk_size = {
            let len = extracted_payloads.len();
            if len < 8 { len } else { len >> 3 }
        };
        let mut compressed_encoded_chunk_vec = Vec::new();
        extracted_payloads
            .par_chunks(chunk_size)
            .map(|chunk| {
                zstd::encode_all(
                    bincode::encode_to_vec(chunk, Self::BINCODE_CONFIG).unwrap().as_slice(),
                    zstd_compression_level,
                )
                .unwrap()
            })
            .collect_into_vec(&mut compressed_encoded_chunk_vec);

        // Put the payloads back into our core
        let _ = std::mem::replace(&mut self.core.payloads, extracted_payloads);

        let header = TTSHeader {
            core_size: compressed_encoded_core.len(),
            payloads_chunk_sizes: compressed_encoded_chunk_vec.iter().map(Vec::len).collect_vec(),
            k_precomputed: self.k_precomputed.clone(),
        };
        let encoded_header = bincode::encode_to_vec(&header, Self::BINCODE_CONFIG)?;

        Ok([
            &Self::MAGIC_NUMBER,
            &(encoded_header.len() as u64).to_le_bytes(),
            encoded_header.as_slice(),
            compressed_encoded_core.as_slice(),
            compressed_encoded_chunk_vec.concat().as_slice(),
        ]
        .concat())
    }

    /// Returns a result containing the tree `core` alongside the `k_precomputed` value
    pub fn from_bytes(data: &[u8]) -> Result<(TaxTreeCore<(), StorePayload>, Vec<usize>), Error> {
        let mut cursor: usize = 0;

        let mut capture_and_advance = |by: usize| {
            let bytes = data.get(cursor..cursor + by).ok_or(Error::MissingBytes(cursor, by));
            cursor += by;
            bytes
        };

        if capture_and_advance(Self::MAGIC_NUMBER.len())? != Self::MAGIC_NUMBER.as_slice() {
            return Err(Error::NotATaxTreeStore);
        }

        let mut header_len: [u8; 8] = [0; 8];
        header_len.copy_from_slice(capture_and_advance(8)?);
        let header_len = u64::from_le_bytes(header_len) as usize;
        let header: TTSHeader =
            bincode::decode_from_slice(capture_and_advance(header_len)?, Self::BINCODE_CONFIG)?.0;

        let mut core: TaxTreeCore<(), StorePayload> = bincode::decode_from_slice(
            zstd::decode_all(capture_and_advance(header.core_size)?)?.as_slice(),
            Self::BINCODE_CONFIG,
        )?
        .0;

        let payloads: Box<[StorePayload]> = header
            .payloads_chunk_sizes
            .into_iter()
            .map(capture_and_advance)
            .collect::<Result<Vec<&[u8]>, Error>>()?
            .into_par_iter()
            .map(|chunk| {
                let decompressed = zstd::decode_all(chunk).unwrap();
                bincode::decode_from_slice::<Vec<StorePayload>, bincode::config::Configuration>(
                    decompressed.as_slice(),
                    Self::BINCODE_CONFIG,
                )
                .unwrap()
                .0
            })
            .flatten_iter()
            .collect::<Vec<StorePayload>>()
            .into_boxed_slice();

        let _ = std::mem::replace(&mut core.payloads, payloads);

        Ok((core, header.k_precomputed))
    }

    pub fn save_to_file(
        &mut self,
        zstd_compression_level: i32,
        multi: Option<&MultiProgress>,
    ) -> Result<(), Error> {
        let spinner = crate::utils::simple_spinner(
            Some(format!(
                "Saving the '{}' gene taxonomic tree to '{}'",
                self.core.gene,
                self.filepath.display()
            )),
            Some(200),
            multi,
        );

        let mut file =
            std::fs::OpenOptions::new().create(true).write(true).truncate(true).open(&self.filepath)?;
        file.write_all(self.to_bytes(zstd_compression_level)?.as_slice())?;
        file.sync_all()?;

        spinner.finish();
        Ok(())
    }

    pub fn load_from_file<Q: AsRef<Path>>(filepath: Q) -> Result<Self, Error> {
        let filepath = filepath.as_ref().to_path_buf();
        let mut file = std::fs::OpenOptions::new().read(true).open(&filepath)?;

        let mut buffer: Vec<u8> = match file.metadata() {
            Ok(metadata) => Vec::with_capacity(metadata.len() as usize),
            Err(e) => {
                eprintln!("Failed to get file metadata, {e}");
                Vec::new()
            }
        };
        let _ = file.read_to_end(&mut buffer)?;

        let (core, k_precomputed) = Self::from_bytes(&buffer)?;

        Ok(Self { core, filepath, k_precomputed })
    }

    /// Very heavy, do this only once if possible
    pub fn load_from_fasta_file<Q: AsRef<Path>>(
        input_filepath: Q,
        store_filepath: Option<Q>,
        gene: impl ToString,
    ) -> Result<Self, Error> {
        let store_filepath = if let Some(fp) = store_filepath {
            fp.as_ref().to_path_buf()
        } else {
            let mut fp = input_filepath.as_ref().to_path_buf();
            let file_name = fp
                .file_name()
                .and_then(|s| s.to_str())
                .map(str::to_owned)
                .unwrap_or_else(|| format!("{}_db", gene.to_string()));
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
        Self::load_from_fasta_string(fasta, store_filepath, gene)
    }

    pub fn load_from_fasta_string(
        fasta: String,
        store_filepath: PathBuf,
        gene: impl ToString,
    ) -> Result<Self, Error> {
        /// Phylogenetic rank name stack (i.e., 'species' is at the top of the stack (end of vector), and 'domain' is at the bottom of the stack (start of vector))
        type Stack = Vec<Box<str>>;
        type StoreNode = Node<()>;

        let mut lines = fasta.lines().enumerate().peekable();

        let highest_rank: OnceCell<Rank> = OnceCell::new();
        let mut leaves_and_stacks: Vec<(StoreNode, Stack)> = Vec::new();
        let mut payloads: Vec<StorePayload> = Vec::new();
        let mut payload_id: usize = 0;

        let spinner = crate::utils::simple_spinner(None::<&str>, Some(200), None);

        let mut string_seq = String::new();

        'outer: while let Some((lidx, text_line)) = lines.next() {
            spinner.set_message(format!("Working on line n°{lidx}"));

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
            let seq: Seq<Iupac> = match Seq::from_str(&string_seq).map_err(|e| Error::BioSeq(e, lidx + 1)) {
                Ok(seq) => seq,
                Err(e) => {
                    eprintln!("{e}");
                    string_seq.clear();
                    continue;
                }
            };

            leaves_and_stacks.push((StoreNode::new_leaf(name, payload_id), stack));
            payloads.push(StorePayload::new(seq, Vec::new()));
            payload_id += 1;
            string_seq.clear();
        }
        spinner.finish();

        // Now we get to the fun part, start at the bottom of the ranks (at species) and move upwards while bundling together
        let spinner =
            crate::utils::simple_spinner(Some("Bundling everything together".to_string()), Some(200), None);

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

        Ok(Self {
            core: TaxTreeCore::new(gene.to_string(), highest_rank, roots, payloads.into_boxed_slice()),
            filepath: store_filepath,
            k_precomputed: Vec::new(),
        })
    }

    pub fn get_filepath(&self) -> &Path {
        &self.filepath
    }

    pub fn set_filepath(
        &mut self,
        filepath: PathBuf,
    ) {
        self.filepath = filepath;
    }

    pub fn shrink(&mut self) {
        self.core.payloads.iter_mut().for_each(|sp| sp.kmer_store_counts_vec.clear());
        self.k_precomputed.clear();
    }

    pub fn compute_and_append_kmer_counts(
        &mut self,
        k: usize,
        nb_resamples: usize,
        multi: Option<&indicatif::MultiProgress>,
    ) {
        let pb = crate::utils::simple_progressbar(
            self.core.payloads.len(),
            format!("computing counts for k={k}"),
            multi,
        );
        self.core.payloads.par_iter_mut().progress_with(pb).for_each_init(
            SmallRng::from_os_rng,
            |rng, sp| {
                let kmer_store_counts =
                    compute_kmer_store_counts_for_fasta_seq(&sp.seq, k, nb_resamples, rng);
                sp.kmer_store_counts_vec.push(kmer_store_counts);
            },
        );
        self.k_precomputed.push(k);
    }

    pub fn compute_and_overwrite_kmer_counts(
        &mut self,
        k: usize,
        nb_resamples: usize,
        multi: Option<&indicatif::MultiProgress>,
    ) {
        let pb = crate::utils::simple_progressbar(
            self.core.payloads.len(),
            format!("computing counts for k={k}"),
            multi,
        );
        self.core.payloads.par_iter_mut().progress_with(pb).for_each_init(
            SmallRng::from_os_rng,
            |rng, sp| {
                let kmer_store_counts =
                    compute_kmer_store_counts_for_fasta_seq(&sp.seq, k, nb_resamples, rng);
                let pre_comp = &mut sp.kmer_store_counts_vec;
                *pre_comp = vec![kmer_store_counts]
            },
        );
        self.k_precomputed = vec![k];
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn to_and_from_bytes_round_trip_test() {
        let leaves = (0..100_usize)
            .map(|payload_id| Node::<()>::new_leaf(Box::from("test-leaf"), payload_id))
            .collect_vec();
        let branch = Node::<()>::new_branch(Box::from("test-branch"), leaves.into_boxed_slice(), ());

        let sp_vec = (0..100_u16)
            .map(|x| StorePayload {
                seq: iupac!("ACTG").to_owned(),
                kmer_store_counts_vec: vec![
                    (
                        unsafe { SparseVec::<u16>::new_unchecked(vec![1, 10, 200], vec![2, 5, 3], 100, x) },
                        0.1
                    );
                    4
                ],
            })
            .collect_vec();

        let core = TaxTreeCore::<(), StorePayload> {
            gene: String::from("test"),
            highest_rank: Rank::Species,
            roots: Box::from([branch]),
            payloads: sp_vec.into_boxed_slice(),
        };

        let mut tree_original = TaxTreeStore {
            core,
            filepath: PathBuf::from("./notimportant.myrtree"),
            k_precomputed: vec![1, 2, 3, 4],
        };

        let tree_reconstructed = {
            let bytes = tree_original.to_bytes(3).unwrap();
            let (core, k_precomputed) = TaxTreeStore::from_bytes(bytes.as_slice()).unwrap();
            TaxTreeStore { core, filepath: PathBuf::from("./notimportant.myrtree"), k_precomputed }
        };

        // Sanity checks
        assert_eq!(tree_original.core.gene, tree_reconstructed.core.gene);
        assert_eq!(tree_original.core.highest_rank, tree_reconstructed.core.highest_rank);
        assert_eq!(tree_original.k_precomputed, tree_reconstructed.k_precomputed);

        assert_eq!(tree_original.core.payloads, tree_reconstructed.core.payloads);
    }
}
