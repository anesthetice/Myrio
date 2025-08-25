use crate::tax::compute_sparse_kmer_counts_for_fasta_seq;
use crate::{data::MyrSeq, tax::clade::Rank};
use crate::{data::SFVec, tax::clade};
use bincode::{Decode, Encode};
use bio_seq::{
    ReverseComplement,
    codec::{dna::Dna, iupac::Iupac},
    error::ParseBioError,
    kmer::Kmer,
    seq::{Seq, SeqSlice},
};
#[cfg(feature = "indicatif")]
use indicatif::{ParallelProgressIterator, ProgressIterator};
use itertools::Itertools;
use myrio_proc::gen_match_k_sparse;
use once_cell::unsync::OnceCell;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};
use std::hint::unreachable_unchecked;
use std::{
    collections::HashMap,
    fs::OpenOptions,
    io::{BufRead, BufReader},
    ops::AddAssign,
    path::Path,
    slice::Iter,
    str::FromStr,
};
use thiserror::Error;

#[derive(Debug, Encode, Decode)]
pub struct TaxTreeStore {
    gene: String,
    roots: Box<[StoreNode]>,
    highest_rank: Rank,
}

#[derive(Debug, Encode, Decode)]
pub enum StoreNode {
    Branch { name: Box<str>, children: Box<[StoreNode]> },
    Leaf { name: Box<str>, seq: Seq<Iupac>, pre_comp: Vec<(usize, SFVec)> },
}

impl TaxTreeStore {
    /// Very heavy, do this only once if possible
    pub fn load_from_fasta_file<Q: AsRef<Path>>(
        filepath: Q,
        gene: impl ToString,
        pre_compute_kcounts: Option<&[usize]>,
        pre_compute_kcounts_max_consecutive_N_before_cutoff: Option<usize>,
    ) -> Result<Self, Error> {
        let file = OpenOptions::new().read(true).open(filepath)?;
        let mut lines = BufReader::new(file).lines().enumerate().peekable();

        let mut highest_rank = Rank::Species;
        let mut leaves_and_stack: Vec<(StoreNode, Vec<Box<str>>)> = Vec::new();
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
        while let Some((lidx, Ok(text_line))) = lines.next() {
            #[cfg(feature = "indicatif")]
            spinner.set_message(format!("working on line n°{lidx}"));

            let (h_rank, _, mut stack) =
                clade::Parsed::from_str(&text_line).map_err(|e| Error::CladeParse(e, lidx + 1))?.uncurl();

            /*
            let expected_stack_size = *stack_size.get_or_init(|| parsed.len() - 1);
            let actual_stack_size = parsed.len() - 1;
            if expected_stack_size != actual_stack_size {
                return Err(Error::RankStackSizeMismatch(expected_stack_size, actual_stack_size, lidx + 1));
            }
            */

            while let Some((_, Ok(next_line))) = lines.peek() {
                if next_line.starts_with('>') {
                    break;
                }
                let dna_line = unsafe { lines.next().unwrap_unchecked().1.unwrap_unchecked() };
                string_seq.push_str(&dna_line);
            }

            let name = unsafe { stack.pop().unwrap_unchecked() }; // safe as an Ok(...) from `clade::parse_str` means the vector isn't empty
            let seq: Seq<Iupac> = Seq::from_str(&string_seq).map_err(|e| Error::BioSeq(e, lidx + 1))?;
            let pre_comp = match pre_compute_kcounts {
                Some(ks) => Vec::with_capacity(ks.len()),
                None => Vec::with_capacity(0),
            };
            leaves_and_stack.push((StoreNode::Leaf { name, seq, pre_comp }, stack));
            string_seq.clear();
        }

        #[cfg(feature = "indicatif")]
        spinner.finish_with_message("Part 1 done");

        if let Some(ks) = pre_compute_kcounts {
            #[allow(non_snake_case)]
            let max_consecutive_N_before_cutoff =
                pre_compute_kcounts_max_consecutive_N_before_cutoff.unwrap_or(4);
            for &k in ks.iter().unique() {
                #[cfg(feature = "indicatif")]
                let bar_count = leaves_and_stack.len() as u64;
                #[cfg(feature = "indicatif")]
                let parallel_iterator = leaves_and_stack.par_iter_mut().progress_count(bar_count);

                #[cfg(not(feature = "indicatif"))]
                let parallel_iterator = leaves_and_stack.par_iter_mut();

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

        Ok(Self { gene: gene.to_string(), roots: Box::from([]), highest_rank })
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
    HighestRankMismatch(Rank, usize, usize),
}
