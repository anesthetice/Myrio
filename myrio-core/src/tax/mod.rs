use std::{
    fs::OpenOptions,
    io::{BufRead, BufReader},
    path::Path,
};

use bio_seq::{codec::iupac::Iupac as DnaIupacCodec, seq::Seq};
use itertools::Itertools;
use once_cell::sync::Lazy;
use regex::Regex;
use thiserror::Error;

use crate::{clustering, clustering::SimFunc, data::MyrSeq};

pub struct TaxTree {
    target: String,
    nodes: Vec<NodeData>,
    children: Vec<Vec<usize>>,
    root: usize,
}

pub enum NodeData {
    Inner { clade: String },
    Tip { name: String, sequence: Seq<DnaIupacCodec> },
}

static DOMAIN_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"a"#).unwrap());

impl TaxTree {
    pub fn load_from_fasta_file<Q: AsRef<Path>>(
        filepath: Q,
        target: String,
    ) -> Result<Self, Error> {
        let mut nodes: Vec<NodeData> = Vec::new();
        let mut children: Vec<Vec<usize>> = Vec::new();

        let insert_entry = |entry: String, sequence: String| {};

        let mut lines = BufReader::new(OpenOptions::new().open(filepath.as_ref())?).lines();

        let mut entry = lines.next().ok_or(Error::NoLines)??;
        if !entry.starts_with(">") {
            return Err(Error::MissingGtSymbol);
        }
        let mut sequence = String::new();

        for line_res in lines {
            let line = line_res?;
            if line.starts_with(">") {}
        }

        todo!()
    }
}

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    IO(#[from] std::io::Error),
    #[error("Invalid or empty FASTA file, contains no lines")]
    NoLines,
    #[error("Invalid FASTA file, first line does not start with `>` symbol")]
    MissingGtSymbol,
}
