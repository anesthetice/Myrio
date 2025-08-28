// Modules
mod func;
pub mod nbin;
pub mod poisson;

// Imports
use itertools::Itertools;
use nbin::{NegativeBinomial as NBinDistr, NegativeBinomialError};
use poisson::{Poisson as PoissonDistr, PoissonError};
use rand::distr::Distribution;

#[derive(Debug, Clone)]
pub enum DiscreteDistribution {
    Cdf { cummul_freqs: &'static [f64] },
    Poisson { inner: PoissonDistr },
    NBin { inner: NBinDistr },
}

impl DiscreteDistribution {
    /// Currently assumes `cummul_freqs` to be valid (i.e., len >= 1 and final value >= 1)
    pub const fn new_cdf(cummul_freqs: &'static [f64]) -> Self {
        Self::Cdf { cummul_freqs }
    }

    pub fn new_poisson(lambda: f64) -> Result<Self, PoissonError> {
        Ok(Self::Poisson { inner: PoissonDistr::new(lambda)? })
    }

    /// The negative binomial distribution is a discrete distribution with two
    /// parameters, `r` and `p`. When `r` is an integer, the negative binomial
    /// distribution can be interpreted as the distribution of the number of
    /// failures in a sequence of Bernoulli trials that continue until `r`
    /// successes occur. `p` is the probability of success in a single Bernoulli trial.
    pub fn new_nbin(
        r: f64,
        p: f64,
    ) -> Result<Self, NegativeBinomialError> {
        Ok(Self::NBin { inner: NBinDistr::new(r, p)? })
    }

    pub fn new_nbin_from_mean_and_std(
        mean: f64,
        std: f64,
    ) -> Result<Self, NegativeBinomialError> {
        let r = mean * mean / (std * std - mean);
        let p = mean / (std * std);
        Ok(Self::NBin { inner: NBinDistr::new(r, p)? })
    }

    pub fn sample_single(
        &self,
        rng: &mut impl rand::Rng,
    ) -> usize {
        match self {
            Self::Cdf { cummul_freqs } => {
                let rand_unif_val: f64 = rng.sample(rand::distr::StandardUniform);
                cummul_freqs.iter().position(|&freq| freq >= rand_unif_val).unwrap()
            }
            Self::Poisson { inner } => inner.sample(rng),
            Self::NBin { inner } => inner.sample(rng),
        }
    }

    pub fn sample_multiple(
        &self,
        n: usize,
        rng: &mut impl rand::Rng,
    ) -> Vec<usize> {
        let mut output: Vec<usize> = Vec::new();
        for _ in 0..n {
            output.push(self.sample_single(rng));
        }
        output
    }
}

pub fn sample_multiple<T>(
    n: usize,
    distr: &impl rand::distr::Distribution<T>,
    rng: &mut impl rand::Rng,
) -> Vec<T> {
    let mut output: Vec<T> = Vec::new();
    for _ in 0..n {
        output.push(rng.sample(distr));
    }
    output
}
