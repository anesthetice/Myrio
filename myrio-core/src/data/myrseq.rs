// Imports
use bio::io::fastq;
use bio_seq::{
    codec::dna::Dna as DnaCodec,
    seq::{Seq, SeqSlice},
};

#[cfg_attr(test, derive(PartialEq))]
#[derive(Clone)]
pub struct MyrSeq {
    /// Sequence identifier
    pub id: String,
    /// Optional description associated with the sequence
    pub description: Option<String>,
    /// Bitpacked sequence of nucleotides
    pub sequence: Seq<DnaCodec>,
    /// Associated 'Phred Quality Score' of each nucleotide
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

    pub fn from_fastq_record_ignore_desc(value: &fastq::Record) -> Self {
        Self {
            id: value.id().to_string(),
            description: None,
            sequence: value.seq().try_into().unwrap(),
            quality: value.qual().iter().map(|a| a.saturating_sub(33)).collect(),
        }
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
}
