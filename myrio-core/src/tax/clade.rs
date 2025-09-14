// Imports
use std::str::FromStr;

use bincode::{Decode, Encode};
use itertools::Itertools;
use once_cell::sync::Lazy;
use regex::Regex;
use thiserror::Error;

// matches `tax={...}`
static TAX_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"tax=\{[^\}]+\}"#).unwrap());

#[derive(Debug, Clone)]
pub struct Parsed {
    pub highest_rank: Rank,
    pub lowest_rank: Rank,
    pub stack: Vec<Box<str>>,
}

impl Parsed {
    pub fn uncurl(self) -> (Rank, Rank, Vec<Box<str>>) {
        (self.highest_rank, self.lowest_rank, self.stack)
    }
}

impl FromStr for Parsed {
    type Err = ParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut view = s;

        match TAX_RE.find(view) {
            Some(_) => view = &view[5..view.len() - 1],
            None => return Err(ParsingError::new("missing `tax={{...}}` delimiter", s)),
        }

        let mut highest_rank: Rank = Rank::Species;
        let mut lowest_rank: Rank = Rank::Domain;

        let stack: Vec<Box<str>> = Rank::all()
            .into_iter()
            .zip_eq(Rank::all_re())
            .filter_map(|(rank, re)| {
                re.find(view).map(|mat| {
                    let name = unsafe { mat.as_str().split_once(":").unwrap_unchecked() }.1.trim();
                    if highest_rank < rank {
                        highest_rank = rank
                    }
                    if lowest_rank > rank {
                        lowest_rank = rank
                    }
                    clean_str(name).into_boxed_str()
                })
            })
            .collect_vec();

        if stack.is_empty() {
            return Err(ParsingError::new("empty entry", s));
        }

        if stack.len() != (highest_rank as usize - lowest_rank as usize + 1) {
            return Err(ParsingError::new(
                format!(
                    "cannot have rank gaps, expected {} elements, got {}",
                    (highest_rank as usize - lowest_rank as usize + 1),
                    stack.len()
                ),
                s,
            ));
        }

        Ok(Self { highest_rank, lowest_rank, stack })
    }
}

#[derive(Debug, Error, PartialEq)]
#[error("Failed to parse string into a list of clade: {msg}; string: '{s}'")]
pub struct ParsingError {
    msg: String,
    s: String,
}

impl ParsingError {
    pub fn new(
        msg: impl ToString,
        s: impl ToString,
    ) -> Self {
        Self { msg: msg.to_string(), s: s.to_string() }
    }
}

#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Encode, Decode)]
#[repr(usize)]
pub enum Rank {
    Domain = 8,
    Kingdom = 7,
    Phylum = 6,
    Class = 5,
    Order = 4,
    Family = 3,
    Genus = 2,
    Species = 1,
}

// matches `domain: Anything really` until it encounters a `,`, `;`, or `}` character, also accepts just `d: Stuff...`
static DOMAIN_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(domain:|d:)[^,;\}]+"#).unwrap());
// et cetera for other ranks
static KINGDOM_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(kingdom:|k:)[^,;\}]+"#).unwrap());
static PHYLUM_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(phylum:|p:)[^,;\}]+"#).unwrap());
static CLASS_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(class:|c:)[^,;\}]+"#).unwrap());
static ORDER_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(order:|o:)[^,;\}]+"#).unwrap());
static FAMILY_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(family:|f:)[^,;\}]+"#).unwrap());
static GENUS_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(genus:|g:)[^,;\}]+"#).unwrap());
static SPECIES_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"(species:|[^s]s:)[^,;\}]+"#).unwrap());

impl Rank {
    pub fn all() -> [Self; 8] {
        [
            Self::Domain,
            Self::Kingdom,
            Self::Phylum,
            Self::Class,
            Self::Order,
            Self::Family,
            Self::Genus,
            Self::Species,
        ]
    }

    pub fn all_re() -> [&'static Lazy<Regex>; 8] {
        [&DOMAIN_RE, &KINGDOM_RE, &PHYLUM_RE, &CLASS_RE, &ORDER_RE, &FAMILY_RE, &GENUS_RE, &SPECIES_RE]
    }

    pub fn above(&self) -> Option<Self> {
        let usize_repr = *self as usize;
        if usize_repr != 8 {
            let rank = unsafe { std::mem::transmute::<usize, Rank>(usize_repr + 1) };
            Some(rank)
        } else {
            None
        }
    }

    pub fn below(&self) -> Option<Self> {
        let usize_repr = *self as usize;
        if usize_repr != 1 {
            let rank = unsafe { std::mem::transmute::<usize, Rank>(usize_repr - 1) };
            Some(rank)
        } else {
            None
        }
    }

    pub fn collect_range_inclusive(
        start: Self,
        stop: Self,
    ) -> Vec<Self> {
        let start_repr_usize = start as usize;
        let stop_repr_usize = stop as usize;

        (start_repr_usize..=stop_repr_usize)
            .map(|i| unsafe { std::mem::transmute::<usize, Rank>(i) })
            .collect_vec()
    }
}

impl std::fmt::Display for Rank {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

fn clean_str(s: &str) -> String {
    let mut output: String = String::with_capacity(s.len());

    let mut capitalize_next_flag: bool = true;

    let mut chars_iter = s.chars().peekable();
    while let Some(chr) = chars_iter.next() {
        if chr.is_whitespace() {
            while let Some(chr) = chars_iter.peek()
                && chr.is_whitespace()
            {
                chars_iter.next();
            }
            output.push(' ');
            capitalize_next_flag = true;
        } else if chr.is_alphabetic() {
            if capitalize_next_flag {
                output.push(chr.to_ascii_uppercase());
                capitalize_next_flag = false;
            } else {
                output.push(chr);
            }
        }
    }
    output
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn clean_str_test() {
        let input: &str = "..gorilla   gorilla Gorilla🫶";
        let expected: &str = "Gorilla Gorilla Gorilla";

        assert_eq!(clean_str(input).as_str(), expected);
    }

    #[test]
    fn parse_test() {
        let input =
            "tax={d:Eukarya, kingdom: Animalia;p: Chordata,, class: The Mammaliàns,   o:A,f:B,g:C,s:D}";
        let expected = vec![
            Box::from("Eukarya"),
            Box::from("Animalia"),
            Box::from("Chordata"),
            Box::from("The Mammaliàns"),
            Box::from("A"),
            Box::from("B"),
            Box::from("C"),
            Box::from("D"),
        ];
        assert_eq!(Parsed::from_str(input).unwrap().stack, expected)
    }
}
