// Imports
use bio::io::fastq;
use bio_seq::{codec::dna::Dna as DnaCodec, seq::Seq};

#[derive(Clone)]
pub struct MyrSeq {
    /// Bitpacked sequence of nucleotides
    pub sequence: Seq<DnaCodec>,
    /// Associated 'Phred Quality Score' of each nucleotide
    pub quality: Vec<u8>,
    /// Optional identifier of the sequence, can be omitted.
    pub id: Option<String>,
    /// Optional description associated with the sequence, can be omitted if it's not required or to save memory.
    pub description: Option<String>,
}

impl From<&fastq::Record> for MyrSeq {
    fn from(value: &fastq::Record) -> Self {
        Self {
            sequence: value.seq().try_into().unwrap(),
            quality: value.qual().iter().map(|a| a.saturating_sub(33)).collect(),
            id: Some(value.id().to_string()),
            description: value.desc().map(str::to_string),
        }
    }
}

impl MyrSeq {
    pub fn new(
        sequence: Seq<DnaCodec>,
        quality: Vec<u8>,
        id: Option<String>,
        description: Option<String>,
    ) -> Self {
        Self { sequence, quality, id, description }
    }

    pub fn from_fastq_record_minimal(value: &fastq::Record) -> Self {
        Self {
            sequence: value.seq().try_into().unwrap(),
            quality: value.qual().iter().map(|a| a.saturating_sub(33)).collect(),
            id: None,
            description: None,
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
        assert_eq!(myrseq.id.as_deref(), Some("1"));
        assert_eq!(myrseq.description.as_deref(), Some("test"));

        // minimal `from`
        let myrseq = MyrSeq::from_fastq_record_minimal(&fastq_record);
        assert_eq!(myrseq.sequence.as_ref(), dna!("ATCG"));
        assert_eq!(myrseq.quality.as_slice(), &[0, 15, 25, 40]);
        assert_eq!(myrseq.id.as_deref(), None);
        assert_eq!(myrseq.description.as_deref(), None);
    }
}
