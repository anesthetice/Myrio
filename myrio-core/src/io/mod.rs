// TODO: add line index to FastqParsingError

// Imports
use std::{
    io::{Read, Write},
    path::Path,
};

use bio_seq::prelude::*;
use itertools::Itertools;
use thiserror::Error;

use crate::MyrSeq;

pub fn read_fastq_from_file<Q: AsRef<Path>>(filepath: Q) -> Result<Vec<MyrSeq>, Error> {
    let filepath: &Path = filepath.as_ref();

    let Some(ext) = filepath.extension().and_then(|s| s.to_str()) else {
        return Err(Error::IO(std::io::Error::new(
            std::io::ErrorKind::InvalidFilename,
            "Failed to retrieve a valid file extension",
        )));
    };

    read_fastq(std::fs::OpenOptions::new().read(true).open(filepath)?, &CompressionMethod::from(ext))
}

pub fn read_fastq<R: Read>(
    input: R,
    compresssion: &CompressionMethod,
) -> Result<Vec<MyrSeq>, Error> {
    let fastq_string = compresssion.decompress_to_string(input)?;

    // Somewhat inspired by the implementation from the `bio` crate, https://docs.rs/bio/latest/src/bio/io/fastq.rs.html#255
    fastq_string
        .lines()
        .tuple_windows::<(_, _, _, _)>()
        .enumerate()
        .step_by(4)
        .map(|(line_idx, (l1, l2, _, l4))| {
            if !l1.starts_with("@") {
                return Err(FastqParsingError::MissingAt(line_idx + 1).into());
            }
            let (id, desc) = match l1[1..].split_once(" ") {
                Some((id, desc)) => (id.to_string(), Some(desc.to_string())),
                None => (l1[1..].to_string(), None),
            };
            let seq: Seq<Dna> =
                Seq::from_str(l2).map_err(|e| Error::Parse(FastqParsingError::BioSeq(e, line_idx + 1)))?;
            let qual: Vec<u8> = l4.as_bytes().iter().map(|val| val.saturating_sub(33)).collect();
            if seq.len() != qual.len() {
                return Err(FastqParsingError::SeqQualLengthMismatch(line_idx + 1).into());
            }
            Ok(MyrSeq::new(id, desc, seq, qual))
        })
        .collect::<Result<Vec<MyrSeq>, Error>>()
}

pub fn write_fastq_to_file<Q: AsRef<Path>>(
    filepath: Q,
    myrseqs: &[MyrSeq],
    compression: &CompressionMethod,
) -> Result<(), Error> {
    let filepath = filepath.as_ref();

    let file = write_fastq(
        std::fs::OpenOptions::new().create(true).write(true).truncate(true).open(filepath)?,
        myrseqs,
        compression,
    )?;
    file.sync_all()?;
    Ok(())
}

pub fn write_fastq<W: Write>(
    output: W,
    myrseqs: &[MyrSeq],
    compression: &CompressionMethod,
) -> Result<W, Error> {
    let data: String = unsafe {
        myrseqs
            .iter()
            .map(|myrseq| {
                [
                    "@",
                    &myrseq.id,
                    if myrseq.description.is_some() { " " } else { "" },
                    myrseq.description.as_deref().unwrap_or(""),
                    "\n",
                    &myrseq.sequence.to_string(),
                    "\n+\n",
                    str::from_utf8_unchecked(
                        &myrseq.quality.iter().map(|q| q.saturating_add(33)).collect_vec(),
                    ),
                ]
                .concat()
            })
            .join("\n")
    };

    compression.compress(data.as_bytes(), output).map_err(Error::from)
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

#[cfg(test)]
mod test {
    use indoc::indoc;

    use super::*;

    #[test]
    fn write_fastq_test() {
        let expected = indoc! {"
            @1
            ACCTTTGGGCCC
            +
            !\"#$%&'()*+,
            @2 test
            ACCTTTGGGC
            +
            !\"#$%&'()*\
        "};

        let input = vec![
            MyrSeq::create("1", None, dna!("ACCTTTGGGCCC"), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
            MyrSeq::create("2", Some("test"), dna!("ACCTTTGGGC"), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        ];

        let output = write_fastq(Vec::new(), &input, &CompressionMethod::None).unwrap();

        assert_eq!(expected.as_bytes(), output.as_slice())
    }

    #[test]
    fn read_fastq_test() {
        let expected = [
            MyrSeq::create("1", None, dna!("ACCTTTGGGCCC"), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
            MyrSeq::create("2", Some("test"), dna!("ACCTTTGGGC"), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        ];

        let input = indoc! {"
            @1
            ACCTTTGGGCCC
            +
            !\"#$%&'()*+,
            @2 test
            ACCTTTGGGC
            +
            !\"#$%&'()*\
        "};

        let output = read_fastq(input.as_bytes(), &CompressionMethod::None).unwrap();

        assert_eq!(expected, output.as_slice())
    }

    #[test]
    fn read_write_fastq_round_trip() {
        let myrseqs = [
            MyrSeq::create("1", None, dna!("ACCTTTGGGCCC"), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
            MyrSeq::create("2", Some("test"), dna!("ACCTTTGGGC"), &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        ];

        let myrseqs_gzip = {
            let compression = CompressionMethod::Gzip(1);
            let data = write_fastq(Vec::new(), &myrseqs, &compression).unwrap();
            read_fastq(data.as_slice(), &compression).unwrap()
        };
        assert_eq!(myrseqs, myrseqs_gzip.as_slice());

        let myrseqs_zstd = {
            let compression = CompressionMethod::Zstd(1);
            let data = write_fastq(Vec::new(), &myrseqs, &compression).unwrap();
            read_fastq(data.as_slice(), &compression).unwrap()
        };
        assert_eq!(myrseqs, myrseqs_zstd.as_slice());
    }
}
