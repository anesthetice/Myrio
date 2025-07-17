// TODO: add line index to FastqParsingError

// Imports
use crate::MyrSeq;
use bio_seq::prelude::*;
use itertools::Itertools;
use std::{
    io::{Read, Write},
    path::Path,
};
use thiserror::Error;

pub fn read_fastq_from_file<Q: AsRef<Path>>(filepath: Q) -> Result<Vec<MyrSeq>, Error> {
    let filepath: &Path = filepath.as_ref();

    let Some(ext) = filepath.extension().and_then(|s| s.to_str()) else {
        return Err(Error::IO(std::io::Error::new(
            std::io::ErrorKind::InvalidFilename,
            "Failed to retrieve a valid file extension",
        )));
    };

    let data = CompressionMethod::from(ext)
        .decompress_to_string(std::fs::OpenOptions::new().read(true).open(filepath)?)?;

    // Somewhat inspired by the implementation from the `bio` crate, https://docs.rs/bio/latest/src/bio/io/fastq.rs.html#255
    let output: Result<Vec<MyrSeq>, Error> = data
        .lines()
        .tuple_windows::<(_, _, _, _)>()
        .enumerate()
        .step_by(4)
        .map(|(line_idx, (l1, l2, _, l4))| {
            if !l1.starts_with("@") {
                return Err(Error::Parse(FastqParsingError::MissingAt(line_idx + 1)));
            }
            let (id, desc) = match l1[1..].split_once(" ") {
                Some((id, desc)) => (id.to_string(), Some(desc.to_string())),
                None => (l1[1..].to_string(), None),
            };
            let seq: Seq<Dna> = Seq::from_str(l2).map_err(|e| FastqParsingError::BioSeq(e, line_idx + 1))?;
            let qual: Vec<u8> = l4.as_bytes().iter().map(|val| val.saturating_sub(33)).collect();
            if seq.len() != qual.len() {
                return Err(Error::Parse(FastqParsingError::SeqQualLengthMismatch(line_idx + 1)));
            }
            Ok(MyrSeq::new(id, desc, seq, qual))
        })
        .collect();

    output
}

pub fn write_fastq_to_file<Q: AsRef<Path>>(
    filepath: Q,
    myrseqs: &[MyrSeq],
    compression: CompressionMethod,
) -> Result<(), Error> {
    let filepath: &Path = filepath.as_ref();

    let data: Vec<u8> = myrseqs
        .iter()
        .map(|myrseq| {
            [
                b"@",
                myrseq.id.as_bytes(),
                b" ",
                myrseq.description.as_deref().map(str::as_bytes).unwrap_or(b""),
                b"\n",
                myrseq.sequence.to_string().as_bytes(),
                b"\n+\n",
                &myrseq.quality.iter().map(|q| q.saturating_add(33)).collect_vec(),
                b"\n",
            ]
            .concat()
        })
        .concat();

    let file_written_to = compression.compress(
        &data,
        std::fs::OpenOptions::new().create(true).write(true).truncate(true).open(filepath)?,
    )?;
    file_written_to.sync_all()?;
    Ok(())
}

pub enum CompressionMethod {
    Gzip(u8),
    Zstd(u8),
    None,
}

impl Default for CompressionMethod {
    fn default() -> Self {
        Self::None
    }
}

impl From<&str> for CompressionMethod {
    fn from(value: &str) -> Self {
        match value.to_ascii_lowercase().as_str() {
            "gz" | "gzip" => Self::Gzip(5),
            "zst" | "zstd" => Self::Zstd(12),
            _ => Self::None,
        }
    }
}

impl CompressionMethod {
    pub fn decompress_to_string<R: Read>(
        &self,
        mut input: R,
    ) -> std::io::Result<String> {
        let mut buffer: String = String::new();
        match self {
            Self::Gzip(_) => {
                let mut decoder = flate2::read::GzDecoder::new(input);
                decoder.read_to_string(&mut buffer)?;
            }
            Self::Zstd(_) => {
                let mut decoder = zstd::Decoder::new(input)?;
                decoder.read_to_string(&mut buffer)?;
            }
            Self::None => {
                input.read_to_string(&mut buffer)?;
            }
        }
        Ok(buffer)
    }

    pub fn compress<W: Write>(
        &self,
        input: &[u8],
        mut output: W,
    ) -> std::io::Result<W> {
        match self {
            Self::Gzip(level) => {
                let mut encoder =
                    flate2::write::GzEncoder::new(output, flate2::Compression::new(*level as u32));
                encoder.write_all(input)?;
                encoder.finish()
            }
            Self::Zstd(level) => {
                let mut encoder = zstd::Encoder::new(output, *level as i32)?;
                encoder.write_all(input)?;
                encoder.finish()
            }
            Self::None => {
                output.write_all(input)?;
                Ok(output)
            }
        }
    }
}

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    IO(#[from] std::io::Error),
    #[error(transparent)]
    Parse(#[from] FastqParsingError),
}

#[derive(Debug, Error)]
pub enum FastqParsingError {
    #[error("Expected an `@` symbol on line {0}")]
    MissingAt(usize),
    #[error(
        "Mismatch between sequence length and number of quality values for the record starting on line {0}"
    )]
    SeqQualLengthMismatch(usize),
    #[error("Bio-seq parsing error for the record starting on line {1}, {0}")]
    BioSeq(ParseBioError, usize),
}
