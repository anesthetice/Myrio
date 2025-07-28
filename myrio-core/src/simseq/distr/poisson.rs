//! Yanked and slightly modfied code from the crate `statrs`, `https://crates.io/crates/statrs`
//!
//! MIT License
//!
//! Copyright (c) 2016 Michael Ma
//!
//! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

use super::func::ln_factorial;
use core::f64;
use thiserror::Error;

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Poisson {
    lambda: f64,
}

#[derive(Copy, Clone, PartialEq, Eq, Error, Debug, Hash)]
pub enum PoissonError {
    #[error("The lambda is NaN, zero, or less than zero")]
    LambdaInvalid,
}

impl Poisson {
    pub fn new(lambda: f64) -> Result<Poisson, PoissonError> {
        if lambda.is_nan() || lambda <= 0.0 {
            Err(PoissonError::LambdaInvalid)
        } else {
            Ok(Poisson { lambda })
        }
    }

    pub fn lambda(&self) -> f64 {
        self.lambda
    }
}

impl rand::distr::Distribution<usize> for Poisson {
    fn sample<R: ::rand::Rng + ?Sized>(
        &self,
        rng: &mut R,
    ) -> usize {
        sample_unchecked(rng, self.lambda) as usize
    }
}

/// Generates one sample from the Poisson distribution either by
/// Knuth's method if lambda < 30.0 or Rejection method PA by
/// A. C. Atkinson from the Journal of the Royal Statistical Society
/// Series C (Applied Statistics) Vol. 28 No. 1. (1979) pp. 29 - 35
/// otherwise
pub(super) fn sample_unchecked<R: rand::Rng + ?Sized>(
    rng: &mut R,
    lambda: f64,
) -> f64 {
    if lambda < 30.0 {
        let limit = (-lambda).exp();
        let mut count = 0.0;
        let mut product: f64 = rng.random();
        while product >= limit {
            count += 1.0;
            product *= rng.random::<f64>();
        }
        count
    } else {
        let c = 0.767 - 3.36 / lambda;
        let beta = f64::consts::PI / (3.0 * lambda).sqrt();
        let alpha = beta * lambda;
        let k = c.ln() - lambda - beta.ln();

        loop {
            let u: f64 = rng.random();
            let x = (alpha - ((1.0 - u) / u).ln()) / beta;
            let n = (x + 0.5).floor();
            if n < 0.0 {
                continue;
            }

            let v: f64 = rng.random();
            let y = alpha - beta * x;
            let temp = 1.0 + y.exp();
            let lhs = y + (v / (temp * temp)).ln();
            let rhs = k + n * lambda.ln() - ln_factorial(n as u64);
            if lhs <= rhs {
                return n;
            }
        }
    }
}
