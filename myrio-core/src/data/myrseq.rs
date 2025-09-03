// Imports
use std::{
    ops::{AddAssign, Range},
    path::Path,
};

use bio_seq::{
    ReverseComplement,
    codec::dna::Dna as DnaCodec,
    seq::{Seq, SeqSlice},
};
use itertools::Itertools;
use myrio_proc::gen_match_k_sparse;
use thiserror::Error;

use crate::{
    constants::Q_TO_BP_CALL_CORRECT_PROB_MAP,
    data::{SFVec, sparse::Float},
};

/// The main data structure used by Myrio, an efficient representation of an FASTQ record
#[cfg_attr(test, derive(PartialEq))]
#[derive(Clone, bincode::Encode, bincode::Decode)]
pub struct MyrSeq {
    /// Sequence identifier
    pub id: String,
    /// Optional description associated with the sequence
    pub description: Option<String>,
    /// Bitpacked sequence of nucleotides
    pub sequence: Seq<DnaCodec>,
    /// Associated 'Quality Score' of each nucleotide
    pub quality: Vec<u8>,
}

impl MyrSeq {
    const BINCODE_CONFIG: bincode::config::Configuration = bincode::config::standard();
    // When `usize` is 64-bit, max{k} = 32 as each nucleotide is repr. by 2 bits
    #[cfg(target_pointer_width = "64")]
    pub const K_SPARSE_VALID_RANGE: Range<usize> = 2..33;
    // When `usize` is 32-bit, max{k} = 16 as each nucleotide is repr. by 2 bits
    #[cfg(target_pointer_width = "32")]
    pub const K_SPARSE_VALID_RANGE: Range<usize> = 2..17;
    #[cfg(target_pointer_width = "64")]
    pub const K_SPARSE_VALID_RANGE_ERROR_MSG: &'static str =
        "For sparse k-mer count maps, only k ∈ {{2, ..., 32}} is currently supported";
    #[cfg(target_pointer_width = "32")]
    pub const K_SPARSE_VALID_RANGE_ERROR_MSG: &'static str =
        "For sparse k-mer count maps, only k ∈ {{2, ..., 16}} is currently supported";

    /// Create a new MyrSeq from owned parameters
    pub fn new(
        id: String,
        desc: Option<String>,
        seq: Seq<DnaCodec>,
        qual: Vec<u8>,
    ) -> Self {
        Self { id, description: desc, sequence: seq, quality: qual }
    }

    /// Create a new MyrSeq from borrowed parameters, useful for creating tests
    pub fn create(
        id: &str,
        desc: Option<&str>,
        seq: &SeqSlice<DnaCodec>,
        qual: &[u8],
    ) -> Self {
        Self {
            id: id.to_string(),
            description: desc.map(str::to_string),
            sequence: seq.to_owned(),
            quality: qual.to_vec(),
        }
    }

    /// Computes the k-mer map. As each k-mer (e.g. `ACTG` if k=4) can be represented as a `usize`, the k-mer uses `usize` keys and returns `f64` values which correspond to the weighted number of a specific k-mer found in the DNA sequence and its reverse complement. Note that only k-mers with a quality score higher than the cutoff are used.
    pub fn compute_kmer_counts(
        &self,
        k: usize,
        cutoff: Float,
    ) -> (SFVec, usize) {
        let conf_score_per_kmer = self
            .quality
            .iter()
            .map(|q| Q_TO_BP_CALL_CORRECT_PROB_MAP[*q as usize])
            .collect_vec()
            .windows(k)
            .map(|vals| vals.iter().product::<Float>())
            .collect_vec();

        let mut nb_hck: usize = 0; // number of high-confidence k-mers
        let mut singles: Vec<usize> = Vec::new();

        macro_rules! body {
            ($seq:expr, $K:expr) => {{
                let nb_kmers = $seq.len() - $K + 1;
                for (idx, (kmer, kmer_rc)) in $seq.kmers::<$K>().zip_eq($seq.to_revcomp().kmers::<$K>()).enumerate() {
                    // Note: `kmer_rc` is not the reverse complement of `kmer`, it's the `idx`-th k-mer of the reverse complement of the sequence.
                    if conf_score_per_kmer[idx] > cutoff {
                        singles.push(usize::from(&kmer));
                        nb_hck += 1;
                    }
                    if conf_score_per_kmer[nb_kmers - 1 - idx] > cutoff {
                        singles.push(usize::from(&kmer_rc));
                        nb_hck += 1;
                    }
                }
            }};
        }
        gen_match_k_sparse!(self.sequence);
        (SFVec::from_unsorted_singles(singles, 1.0, 4_usize.pow(k as u32), 0.0), nb_hck)
    }

    pub fn pre_process(
        input: Vec<MyrSeq>,
        min_mean_qual: Float,
        max_qual: usize,
        min_length: usize,
    ) -> Vec<MyrSeq> {
        todo!()
    }

    pub fn encode_vec<W: std::io::Write>(
        myrseqs: &[MyrSeq],
        compression_level: i32,
        output: W,
    ) -> Result<W, Error> {
        let mut encoder = zstd::Encoder::new(output, compression_level)?;
        bincode::encode_into_std_write(myrseqs, &mut encoder, Self::BINCODE_CONFIG)?;
        encoder.finish().map_err(Error::from)
    }

    pub fn encode_vec_to_file<Q: AsRef<Path>>(
        filepath: Q,
        myrseqs: &[MyrSeq],
        compression_level: i32,
    ) -> Result<(), Error> {
        let file = Self::encode_vec(
            myrseqs,
            compression_level,
            std::fs::OpenOptions::new().create(true).write(true).truncate(true).open(filepath)?,
        )?;

        file.sync_all().map_err(Error::from)
    }

    pub fn decode_vec<R: std::io::Read>(input: R) -> Result<Vec<MyrSeq>, Error> {
        let mut decoder = zstd::Decoder::new(input)?;
        bincode::decode_from_std_read(&mut decoder, Self::BINCODE_CONFIG).map_err(Error::from)
    }

    pub fn decode_vec_from_file<Q: AsRef<Path>>(filepath: Q) -> Result<Vec<MyrSeq>, Error> {
        Self::decode_vec(std::fs::OpenOptions::new().read(true).open(filepath)?)
    }
}

impl core::fmt::Debug for MyrSeq {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        if f.alternate() {
            write!(
                f,
                "MyrSeq {{\n    sequence: {}\n    quality: {:?}\n    id: {:?}\n    description: {:?}\n}}",
                self.sequence, self.quality, self.id, self.description
            )
        } else {
            write!(
                f,
                "MyrSeq {{seq:{}, qual:{:?}, id:{:?}, desc:{:?}}}",
                self.sequence, self.quality, self.id, self.description
            )
        }
    }
}

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    BincodeDecode(#[from] bincode::error::DecodeError),
    #[error(transparent)]
    BincodeEncode(#[from] bincode::error::EncodeError),
    #[error(transparent)]
    IO(#[from] std::io::Error),
}

#[cfg(test)]
mod test {
    use bio_seq::prelude::*;

    use super::MyrSeq;

    #[test]
    pub fn encode_decode_round_trip_test() {
        let myrseqs = [
            MyrSeq::create("1", None, dna!("ACCTTTGGGCCC"), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
            MyrSeq::create("2", Some("test"), dna!("ACCTTTGGGC"), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        ];

        let reconstructed_myrseqs = {
            let data = MyrSeq::encode_vec(&myrseqs, 3, Vec::new()).unwrap();
            MyrSeq::decode_vec(data.as_slice()).unwrap()
        };

        assert_eq!(myrseqs, reconstructed_myrseqs.as_slice())
    }
}
