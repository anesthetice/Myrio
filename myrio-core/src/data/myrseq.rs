// Imports
use std::{
    collections::HashMap,
    ops::{AddAssign, Range},
    path::Path,
};

use bio::io::fastq;
use bio_seq::{
    ReverseComplement,
    codec::dna::Dna as DnaCodec,
    seq::{Seq, SeqSlice},
};
use itertools::Itertools;
use myrio_proc::{gen_match_k_dense, gen_match_k_sparse};
use ndarray::Array1;
use thiserror::Error;

use super::sparse::SparseFloatVec;
use crate::constants::Q_TO_BP_CALL_CORRECT_PROB_MAP;

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

impl From<&fastq::Record> for MyrSeq {
    fn from(value: &fastq::Record) -> Self {
        Self {
            id: value.id().to_string(),
            description: value.desc().map(str::to_string),
            sequence: value.seq().try_into().unwrap(),
            quality: value.qual().iter().map(|a| a.saturating_sub(33)).collect(),
        }
    }
}

impl MyrSeq {
    pub const K_DENSE_VALID_RANGE: Range<usize> = 2..10;
    pub const K_DENSE_VALID_RANGE_ERROR_MSG: &'static str =
        "For dense k-mer count maps, only k ∈ {{2, ..., 9}} is currently supported";
    pub const K_SPARSE_VALID_RANGE: Range<usize> = 2..43;
    pub const K_SPARSE_VALID_RANGE_ERROR_MSG: &'static str =
        "For sparse k-mer count maps, only k ∈ {{2, ..., 42}} is currently supported";

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

    pub fn compute_dense_kmer_counts(
        &self,
        k: usize,
        cutoff: f64,
    ) -> Result<(Array1<f64>, usize), Error> {
        let conf_score_per_kmer = self
            .quality
            .iter()
            .map(|q| Q_TO_BP_CALL_CORRECT_PROB_MAP[*q as usize])
            .collect_vec()
            .windows(k)
            .map(|vals| vals.iter().product::<f64>())
            .collect_vec();

        macro_rules! body {
            ($seq:expr, $K:expr) => {{
                let nb_kmers = $seq.len() - $K + 1;
                let mut nb_hck: usize = 0; // number of high-quality k-mers
                let mut map = Array1::<f64>::zeros(4_usize.pow(k as u32));

                for (idx, (kmer, kmer_rc)) in $seq.kmers::<$K>().zip_eq($seq.to_revcomp().kmers::<$K>()).enumerate() {
                    // Note: `kmer_rc` is not the reverse complement of `kmer`, it's the `idx`-th k-mer of the reverse complement of the sequence; cleaner implementation if `KmerIter` supported `.rev()` method but it doesn't unfortunately.
                    if conf_score_per_kmer[idx] > cutoff {
                        map[usize::from(&kmer)] += 1.0;
                        nb_hck += 1;
                    }
                    if conf_score_per_kmer[nb_kmers - 1 - idx] > cutoff {
                        map[usize::from(&kmer_rc)] += 1.0;
                        nb_hck += 1;
                    }
                }
                Ok((map, nb_hck))
            }};
        }
        gen_match_k_dense!(self.sequence)
    }

    /// Computes the k-mer map. As each k-mer (e.g. `ACTG` if k=4) can be represented as a `usize`, the k-mer uses `usize` keys and returns `f64` values which correspond to the weighted number of a specific k-mer found in the DNA sequence and its reverse complement. Note that only k-mers with a quality score higher than the cutoff are used.
    pub fn compute_sparse_kmer_counts(
        &self,
        k: usize,
        cutoff: f64,
    ) -> Result<(SparseFloatVec, usize), Error> {
        let conf_score_per_kmer = self
            .quality
            .iter()
            .map(|q| Q_TO_BP_CALL_CORRECT_PROB_MAP[*q as usize])
            .collect_vec()
            .windows(k)
            .map(|vals| vals.iter().product::<f64>())
            .collect_vec();

        macro_rules! body {
            ($seq:expr, $K:expr) => {{
                let nb_kmers = $seq.len() - $K + 1;
                let mut nb_hck: usize = 0; // number of high-quality k-mers
                let mut map = SparseFloatVec::default();

                for (idx, (kmer, kmer_rc)) in $seq.kmers::<$K>().zip_eq($seq.to_revcomp().kmers::<$K>()).enumerate() {
                    // Note: `kmer_rc` is not the reverse complement of `kmer`, it's the `idx`-th k-mer of the reverse complement of the sequence; cleaner implementation if `KmerIter` supported `.rev()` method but it doesn't unfortunately.
                    if conf_score_per_kmer[idx] > cutoff {
                        map.get_or_insert_mut(usize::from(&kmer)).add_assign(1.0);
                        nb_hck += 1;

                    }
                    if conf_score_per_kmer[nb_kmers - 1 - idx] > cutoff {
                        map.get_or_insert_mut(usize::from(&kmer_rc)).add_assign(1.0);
                        nb_hck += 1;
                    }
                }
                Ok((map, nb_hck))
            }};
        }
        gen_match_k_sparse!(self.sequence)
    }

    pub fn from_fastq_record_ignore_desc(value: &fastq::Record) -> Self {
        Self {
            id: value.id().to_string(),
            description: None,
            sequence: value.seq().try_into().unwrap(),
            quality: value.qual().iter().map(|a| a.saturating_sub(33)).collect(),
        }
    }

    pub fn encode_vec<W: std::io::Write>(
        myrseqs: &[MyrSeq],
        compression_level: i32,
        output: W,
    ) -> Result<W, Error> {
        let config = bincode::config::standard();
        let mut encoder = zstd::Encoder::new(output, compression_level)?;
        bincode::encode_into_std_write(myrseqs, &mut encoder, config)?;
        encoder.finish().map_err(Error::from)
    }

    pub fn encode_vec_to_file<Q: AsRef<Path>>(
        filepath: Q,
        myrseqs: &[MyrSeq],
        compression_level: i32,
    ) -> Result<(), Error> {
        let filepath = filepath.as_ref();

        let file = Self::encode_vec(
            myrseqs,
            compression_level,
            std::fs::OpenOptions::new().create(true).write(true).truncate(true).open(filepath)?,
        )?;

        file.sync_all().map_err(Error::from)
    }

    pub fn decode_vec<R: std::io::Read>(input: R) -> Result<Vec<MyrSeq>, Error> {
        let config = bincode::config::standard();
        let mut decoder = zstd::Decoder::new(input)?;
        bincode::decode_from_std_read(&mut decoder, config).map_err(Error::from)
    }

    pub fn decode_vec_from_file<Q: AsRef<Path>>(filepath: Q) -> Result<Vec<MyrSeq>, Error> {
        let filepath: &Path = filepath.as_ref();
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
    BincodeDecodeError(#[from] bincode::error::DecodeError),
    #[error(transparent)]
    BincodeEncodeError(#[from] bincode::error::EncodeError),
    #[error(transparent)]
    IO(#[from] std::io::Error),
    #[error("{0}")]
    InvalidKmerSize(&'static str),
}

#[cfg(test)]
mod test {
    use bio::io::fastq;
    use bio_seq::prelude::*;

    use super::MyrSeq;

    #[test]
    pub fn test_myrseq_from_fastq_record() {
        let fastq_record = fastq::Record::with_attrs("1", Some("test"), b"ATCG", b"!0:I");

        // default `from`
        let myrseq = MyrSeq::from(&fastq_record);
        assert_eq!(myrseq.sequence.as_ref(), dna!("ATCG"));
        assert_eq!(myrseq.quality.as_slice(), &[0, 15, 25, 40]);
        assert_eq!(myrseq.id.as_str(), "1");
        assert_eq!(myrseq.description.as_deref(), Some("test"));

        // minimal `from`
        let myrseq = MyrSeq::from_fastq_record_ignore_desc(&fastq_record);
        assert_eq!(myrseq.sequence.as_ref(), dna!("ATCG"));
        assert_eq!(myrseq.quality.as_slice(), &[0, 15, 25, 40]);
        assert_eq!(myrseq.id.as_str(), "1");
        assert_eq!(myrseq.description, None);
    }

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
