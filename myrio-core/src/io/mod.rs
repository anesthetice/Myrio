// TODO: add line index to FastqParsingError

// Imports
use crate::MyrSeq;
use bio_seq::prelude::*;
use itertools::Itertools;
use std::{io::Read, path::Path};
use thiserror::Error;

pub fn read_fastq<Q: AsRef<Path>>(filepath: Q) -> Result<Vec<MyrSeq>, Error> {
    let filepath: &Path = filepath.as_ref();

    let Some(ext) = filepath.extension().and_then(|s| s.to_str()) else {
        return Err(Error::InvalidFilepath("missing valid extension"));
    };

    let mut file = std::fs::OpenOptions::new().read(true).open(filepath)?;
    let mut data: String = String::new();
    match ext {
        "fastq" => file.read_to_string(&mut data)?,
        "zst" => {
            let mut decoder = zstd::Decoder::new(file)?;
            decoder.read_to_string(&mut data)?
        }
        "gz" => {
            let mut decoder = flate2::read::GzDecoder::new(file);
            decoder.read_to_string(&mut data)?
        }
        _ => return Err(Error::InvalidFilepath("unknown extension")),
    };

    // Somewhat inspired by the implementation from the `bio` crate, https://docs.rs/bio/latest/src/bio/io/fastq.rs.html#255
    let output: Result<Vec<MyrSeq>, Error> = data
        .lines()
        .tuple_windows::<(_, _, _, _)>()
        .step_by(4)
        .map(|(l1, l2, _, l4)| {
            if !l1.starts_with("@") {
                return Err(Error::Parse(FastqParsingError::MissingAt));
            }
            let (id, desc) = match l1[1..].split_once(" ") {
                Some((id, desc)) => (id.to_string(), Some(desc.to_string())),
                None => (l1[1..].to_string(), None),
            };
            let seq: Seq<Dna> = Seq::from_str(l2).map_err(FastqParsingError::BioSeq)?;
            let qual: Vec<u8> = l4.as_bytes().iter().map(|val| val.saturating_sub(33)).collect();
            if seq.len() != qual.len() {
                return Err(Error::Parse(FastqParsingError::Mismatch));
            }
            Ok(MyrSeq::new(id, desc, seq, qual))
        })
        .collect();

    output
}

#[derive(Debug, Error)]
pub enum Error {
    #[error("Invalid filepath: {0}")]
    InvalidFilepath(&'static str),
    #[error(transparent)]
    IO(#[from] std::io::Error),
    #[error(transparent)]
    Parse(#[from] FastqParsingError),
}

#[derive(Debug, Error)]
pub enum FastqParsingError {
    #[error("Missing `@` symbol")]
    MissingAt,
    #[error("Mismatch between sequence length and number of qualities")]
    Mismatch,
    #[error(transparent)]
    BioSeq(#[from] ParseBioError),
}
