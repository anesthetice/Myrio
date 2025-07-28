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

#![allow(clippy::excessive_precision)]

use core::f64;

pub(super) fn ln_factorial(x: u64) -> f64 {
    let x = x as usize;
    FCACHE.get(x).map_or_else(|| ln_gamma(x as f64 + 1.0), |&fac| fac.ln())
}

/// The maximum factorial representable
/// by a 64-bit floating point without
/// overflowing
pub(super) const MAX_FACTORIAL: usize = 170;

// Initialization for pre-computed cache of 171 factorial
// values 0!...170!
const FCACHE: [f64; MAX_FACTORIAL + 1] = {
    let mut fcache = [1.0; MAX_FACTORIAL + 1];

    // `const` only allow while loops
    let mut i = 1;
    while i < MAX_FACTORIAL + 1 {
        fcache[i] = fcache[i - 1] * i as f64;
        i += 1;
    }

    fcache
};

/// Computes the logarithm of the gamma function
/// with an accuracy of 16 floating point digits.
/// The implementation is derived from
/// "An Analysis of the Lanczos Gamma Approximation",
/// Glendon Ralph Pugh, 2004 p. 116
pub(super) fn ln_gamma(x: f64) -> f64 {
    if x < 0.5 {
        let s = GAMMA_DK.iter().enumerate().skip(1).fold(GAMMA_DK[0], |s, t| s + t.1 / (t.0 as f64 - x));

        LN_PI
            - (f64::consts::PI * x).sin().ln()
            - s.ln()
            - LN_2_SQRT_E_OVER_PI
            - (0.5 - x) * ((0.5 - x + GAMMA_R) / f64::consts::E).ln()
    } else {
        let s =
            GAMMA_DK.iter().enumerate().skip(1).fold(GAMMA_DK[0], |s, t| s + t.1 / (x + t.0 as f64 - 1.0));

        s.ln() + LN_2_SQRT_E_OVER_PI + (x - 0.5) * ((x - 0.5 + GAMMA_R) / f64::consts::E).ln()
    }
}

/// Auxiliary variable when evaluating the `gamma_ln` function
const GAMMA_R: f64 = 10.900511;

/// Polynomial coefficients for approximating the `gamma_ln` function
const GAMMA_DK: &[f64] = &[
    2.48574089138753565546e-5,
    1.05142378581721974210,
    -3.45687097222016235469,
    4.51227709466894823700,
    -2.98285225323576655721,
    1.05639711577126713077,
    -1.95428773191645869583e-1,
    1.70970543404441224307e-2,
    -5.71926117404305781283e-4,
    4.63399473359905636708e-6,
    -2.71994908488607703910e-9,
];

/// Constant value for `ln(pi)`
pub(super) const LN_PI: f64 = 1.1447298858494001741434273513530587116472948129153;

/// Constant value for `ln(2 * sqrt(e / pi))`
pub(super) const LN_2_SQRT_E_OVER_PI: f64 = 0.6207822376352452223455184457816472122518527279025978;
