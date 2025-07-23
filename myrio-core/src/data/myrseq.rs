// Imports
use std::{
    collections::HashMap,
    ops::{AddAssign, Range},
    path::Path,
};

use bio::io::fastq;
use bio_seq::{
    codec::dna::Dna as DnaCodec,
    prelude::*,
    seq::{Seq, SeqSlice},
};
use itertools::Itertools;
use thiserror::Error;

use crate::constants::Q_TO_BP_CALL_CORRECT_PROB_MAP;

/// The main data structure used by Myrio
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
    pub const K_VALID_RANGE: Range<usize> = 2..43;
    pub const K_VALID_RANGE_ERROR_MSG: &'static str = "Only k ∈ {{2, ..., 42}} is currently supported";

    pub fn new(
        id: String,
        desc: Option<String>,
        seq: Seq<DnaCodec>,
        qual: Vec<u8>,
    ) -> Self {
        Self { id, description: desc, sequence: seq, quality: qual }
    }

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

    pub fn get_kmer_map(
        &self,
        k: usize,
        cutoff: f64,
    ) -> Result<(HashMap<usize, f64>, f64), Error> {
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

                #[allow(non_snake_case)]
                let mut nb_HCS_weighted: f64 = 0.0;
                let mut map: HashMap<usize, f64> = HashMap::new();

                for (idx, (kmer, kmer_rc)) in $seq.kmers::<$K>().zip_eq($seq.to_revcomp().kmers::<$K>()).enumerate() {
                    // Note: `kmer_rc` is not the reverse complement of `kmer`, it's the `idx`-th k-mer of the reverse complement of the sequence; cleaner implementation if `KmerIter` supported `.rev()` method but it doesn't unfortunately.
                    let conf_score = conf_score_per_kmer[idx];
                    if conf_score > cutoff {
                        let val_ref = map.entry(usize::from(&kmer)).or_default();
                        nb_HCS_weighted += conf_score;
                        val_ref.add_assign(conf_score)
                    }
                    let conf_score = conf_score_per_kmer[nb_kmers - 1 - idx];
                    if conf_score > cutoff {
                        let val_ref = map.entry(usize::from(&kmer_rc)).or_default();
                        nb_HCS_weighted += conf_score;
                        val_ref.add_assign(conf_score)
                    }
                }
                Ok((map, nb_HCS_weighted))
            }};
        }

        match k {
            2 => body!(self.sequence, 2),
            3 => body!(self.sequence, 3),
            4 => body!(self.sequence, 4),
            5 => body!(self.sequence, 5),
            6 => body!(self.sequence, 6),
            7 => body!(self.sequence, 7),
            8 => body!(self.sequence, 8),
            9 => body!(self.sequence, 9),
            10 => body!(self.sequence, 10),
            11 => body!(self.sequence, 11),
            12 => body!(self.sequence, 12),
            13 => body!(self.sequence, 13),
            14 => body!(self.sequence, 14),
            15 => body!(self.sequence, 15),
            16 => body!(self.sequence, 16),
            17 => body!(self.sequence, 17),
            18 => body!(self.sequence, 18),
            19 => body!(self.sequence, 19),
            20 => body!(self.sequence, 20),
            21 => body!(self.sequence, 21),
            22 => body!(self.sequence, 22),
            23 => body!(self.sequence, 23),
            24 => body!(self.sequence, 24),
            25 => body!(self.sequence, 25),
            26 => body!(self.sequence, 26),
            27 => body!(self.sequence, 27),
            28 => body!(self.sequence, 28),
            29 => body!(self.sequence, 29),
            30 => body!(self.sequence, 30),
            31 => body!(self.sequence, 31),
            32 => body!(self.sequence, 32),
            33 => body!(self.sequence, 33),
            34 => body!(self.sequence, 34),
            35 => body!(self.sequence, 35),
            36 => body!(self.sequence, 36),
            37 => body!(self.sequence, 37),
            38 => body!(self.sequence, 38),
            39 => body!(self.sequence, 39),
            40 => body!(self.sequence, 40),
            41 => body!(self.sequence, 41),
            42 => body!(self.sequence, 42),
            _ => Err(Error::InvalidKmerSize),
        }
    }

    pub fn get_kmer_map_or_panic(
        &self,
        k: usize,
        cutoff: f64,
    ) -> (HashMap<usize, f64>, f64) {
        self.get_kmer_map(k, cutoff).expect(Self::K_VALID_RANGE_ERROR_MSG)
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
    #[error("{}", MyrSeq::K_VALID_RANGE_ERROR_MSG)]
    InvalidKmerSize,
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
