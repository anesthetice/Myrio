// Imports

use itertools::Itertools;
use once_cell::sync::Lazy;
use regex::Regex;
use thiserror::Error;

// matches `tax={...}`
static TAX_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r#"tax=\{[^\}]+\}"#).unwrap());

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

#[derive(Debug, Clone, PartialEq)]
pub struct Clade {
    pub name: String,
    pub rank: Rank,
}

impl Clade {
    pub fn new(
        name: &str,
        rank: Rank,
    ) -> Self {
        Self { name: name.to_string(), rank }
    }

    pub fn parse_str(s: &str) -> Result<Vec<Self>, ParsingError> {
        let mut view = s;

        match TAX_RE.find(view) {
            Some(_) => view = &view[5..view.len() - 1],
            None => return Err(ParsingError::new("missing `tax={{...}}` delimiter", s)),
        }

        let output =
            [&DOMAIN_RE, &KINGDOM_RE, &PHYLUM_RE, &CLASS_RE, &ORDER_RE, &FAMILY_RE, &GENUS_RE, &SPECIES_RE]
                .into_iter()
                .zip_eq([
                    Rank::Domain,
                    Rank::Kingdom,
                    Rank::Phylum,
                    Rank::Class,
                    Rank::Order,
                    Rank::Family,
                    Rank::Genus,
                    Rank::Species,
                ])
                .filter_map(|(re, rank)| {
                    re.find(view).map(|mat| {
                        let name = unsafe { mat.as_str().split_once(":").unwrap_unchecked() }.1.trim();
                        Self::new(name, rank)
                    })
                })
                .collect_vec();

        if output.is_empty() {
            return Err(ParsingError::new("empty entry", s));
        }

        if output.len() != output[0].rank as usize {
            return Err(ParsingError::new(
                format!(
                    "cannot have rank gaps, expected {} elements due to `{:?}` being set, got only {} instead",
                    output[0].rank as usize,
                    output[0].rank,
                    output.len()
                ),
                s,
            ));
        }

        Ok(output)
    }
}

#[derive(Debug, Error, PartialEq)]
#[error("Failed to parse string into a list of clade: {msg}\nstring: '{s}'")]
pub struct ParsingError<'a> {
    msg: String,
    s: &'a str,
}

impl<'a> ParsingError<'a> {
    pub fn new(
        msg: impl ToString,
        s: &'a str,
    ) -> Self {
        Self { msg: msg.to_string(), s }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn parse_test() {
        let input =
            "tax={d:Eukarya, kingdom: Animalia;p: Chordata,, class: The Mammaliàns,   o:A,f:B,g:C,s:D}";
        let expected = vec![
            Clade::new("Eukarya", Rank::Domain),
            Clade::new("Animalia", Rank::Kingdom),
            Clade::new("Chordata", Rank::Phylum),
            Clade::new("The Mammaliàns", Rank::Class),
            Clade::new("A", Rank::Order),
            Clade::new("B", Rank::Family),
            Clade::new("C", Rank::Genus),
            Clade::new("D", Rank::Species),
        ];
        assert_eq!(Clade::parse_str(input).unwrap(), expected)
    }
}
