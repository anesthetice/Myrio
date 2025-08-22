use std::{
    fs::OpenOptions,
    io::{BufRead, BufReader},
    path::Path,
};

use bincode::{Decode, Encode};
use thiserror::Error;

use crate::data::SFVec;

// todo: use smallvec instead?
#[derive(Debug, Encode, Decode)]
pub struct TaxTreeStore {
    gene: String,
    domains: Vec<(String, Vec<usize>)>,
    kingdoms: Vec<(String, Vec<usize>)>,
    phylums: Vec<(String, Vec<usize>)>,
    classes: Vec<(String, Vec<usize>)>,
    orders: Vec<(String, Vec<usize>)>,
    families: Vec<(String, Vec<usize>)>,
    genuses: Vec<(String, Vec<usize>)>,
    species: Vec<(String, SFVec)>,
}

impl TaxTreeStore {
    /// Very heavy, do this only once if possible
    pub fn load_from_fasta_file<Q: AsRef<Path>>(
        filepath: Q,
        target: String,
    ) -> Result<Self, Error> {
        let file = OpenOptions::new().read(true).open(filepath)?;
        let mut lines = BufReader::new(file).lines();

        while let Some(text_line) = lines.next() {
            let mut full_sequence = String::new();
            crate::tax::clade::Clade::parse_str(text_line?.as_str()).unwrap();
            //while let Some(dna_line)
        }

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
