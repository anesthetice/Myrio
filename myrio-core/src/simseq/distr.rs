use itertools::Itertools;
use rand_distr::{Distribution, Normal as NormalRD, Poisson as PoissonRD};

#[derive(Debug, Clone)]
pub enum UsizeDistribution {
    Cdf { cummul_freqs: &'static [f64] },
    Poisson { inner: PoissonRD<f64> },
}

impl UsizeDistribution {
    /// Currently assumes `cummul_freqs` to be valid (i.e., len >= 1 and final value >= 1)
    pub fn new_cdf(cummul_freqs: &'static [f64]) -> Self {
        Self::Cdf { cummul_freqs }
    }

    pub fn new_poisson(lambda: f64) -> Result<Self, rand_distr::PoissonError> {
        Ok(Self::Poisson { inner: PoissonRD::new(lambda)? })
    }

    pub fn sample_single(
        &self,
        rng: &mut impl rand::Rng,
    ) -> usize {
        match self {
            Self::Cdf { cummul_freqs } => {
                let rand_unif_val: f64 = rng.sample(rand::distr::StandardUniform);
                cummul_freqs.iter().find_position(|freq| **freq >= rand_unif_val).unwrap().0
            }
            Self::Poisson { inner } => inner.sample(rng).round() as usize,
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

#[derive(Debug, Clone)]
pub enum FloatDistribution {
    Normal { inner: NormalRD<f64> },
}

impl FloatDistribution {
    pub fn new_normal(
        mean: f64,
        std_dev: f64,
    ) -> Result<Self, rand_distr::NormalError> {
        Ok(Self::Normal { inner: NormalRD::new(mean, std_dev)? })
    }

    pub fn sample_single(
        &self,
        rng: &mut impl rand::Rng,
    ) -> f64 {
        match self {
            Self::Normal { inner } => inner.sample(rng),
        }
    }

    pub fn sample_multiple(
        &self,
        n: usize,
        rng: &mut impl rand::Rng,
    ) -> Vec<f64> {
        let mut output: Vec<f64> = Vec::new();
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
