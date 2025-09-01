// Modules
pub mod constants;
pub mod distr;

// Imports
use std::ops::SubAssign;

use bio_seq::prelude::*;
use distr::{DiscreteDistribution, sample_multiple};
use itertools::Itertools;
use myrio_core::data::MyrSeq;
use rand::{distr::Distribution, seq::IndexedRandom};
use thiserror::Error;

use crate::simseq::constants::{
    MAX_Q_SCORE, MIN_Q_SCORE, OBSERVED_Q_SCORE_BLOCK_SIZE_CUMMUL_FREQ, OBSERVED_Q_SCORE_CUMMUL_FREQ,
};

#[derive(Debug)]
pub struct Generator {
    pub coding_to_template_ratio_bounds: (f64, f64),
    pub q_score_distr: DiscreteDistribution,
    pub q_score_block_size_distr: DiscreteDistribution,
    pub indel_insertion_snp_error_weights: (f64, f64, f64),
    pub min_window_size: usize,
}

impl Default for Generator {
    fn default() -> Self {
        Self {
            coding_to_template_ratio_bounds: (0.35, 0.65),
            q_score_distr: DiscreteDistribution::new_cdf(&OBSERVED_Q_SCORE_CUMMUL_FREQ),
            q_score_block_size_distr: DiscreteDistribution::new_cdf(&OBSERVED_Q_SCORE_BLOCK_SIZE_CUMMUL_FREQ),
            indel_insertion_snp_error_weights: (0.3, 0.3, 0.4),
            min_window_size: 150,
        }
    }
}

impl Generator {
    /// `x = 0` is at index `3`
    const ADJUSTMENT: usize = 3;
    /// Represents a simple distribution of how q-scores are offset within a block of q-scores from the main q-score
    const Q_SCORE_INNER_BLOCK_OFFSET_DISTR: DiscreteDistribution =
        DiscreteDistribution::new_cdf(&[0.05, 0.15, 0.30, 0.70, 0.85, 0.95, 1.0]);

    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_coding_to_template_ratio_bounds(
        self,
        bounds: (f64, f64),
    ) -> Result<Self, Error> {
        if bounds.0 > bounds.1 || bounds.0 < 0.0 || bounds.0 > 1.0 || bounds.1 < 0.0 || bounds.1 > 1.0 {
            Err(Error::InvalidCTRBounds)
        } else {
            Ok(Self { coding_to_template_ratio_bounds: bounds, ..self })
        }
    }

    pub fn with_q_score_distr(
        self,
        distr: DiscreteDistribution,
    ) -> Self {
        Self { q_score_distr: distr, ..self }
    }

    pub fn with_q_score_block_size_distr(
        self,
        distr: DiscreteDistribution,
    ) -> Self {
        Self { q_score_block_size_distr: distr, ..self }
    }

    pub fn with_indel_insertion_snp_error_weights(
        self,
        weights: (f64, f64, f64),
    ) -> Result<Self, Error> {
        if weights.0.is_sign_negative()
            || !weights.0.is_finite()
            || weights.1.is_sign_negative()
            || !weights.1.is_finite()
            || weights.2.is_sign_negative()
            || !weights.2.is_finite()
        {
            Err(Error::InvalidSequencingErrorWeights)
        } else {
            Ok(Self { indel_insertion_snp_error_weights: weights, ..self })
        }
    }

    pub fn with_min_window_size(
        self,
        size: usize,
    ) -> Self {
        Self { min_window_size: size, ..self }
    }

    pub fn generate_core_sequence(
        length: usize,
        rng: &mut impl rand::Rng,
    ) -> Seq<Dna> {
        let mut buffer = vec![0_u8; ((usize::BITS as usize) >> 3) * length];
        rng.fill_bytes(&mut buffer);
        Seq::from_raw(length, bytemuck::cast_slice(&buffer)).unwrap()
    }

    pub fn generate_pseudo_amplicon(
        &self,
        length: usize,
        amount: usize,
        id: &str,
        rng: &mut impl rand::Rng,
    ) -> Vec<MyrSeq> {
        // We generate our core DNA forward amplicon sequence, and its reverse complement
        let core_seq: Seq<Dna> = Self::generate_core_sequence(length, rng);
        let core_seq_rc = core_seq.to_revcomp();

        // We generate the ratio of sequences corresponding to the coding strand and the template (reverse complement) strand
        let forward_seq_ratio = rand::distr::uniform::Uniform::new(
            self.coding_to_template_ratio_bounds.0,
            self.coding_to_template_ratio_bounds.1,
        )
        .unwrap()
        .sample(&mut *rng);

        let w_size_distr =
            DiscreteDistribution::new_nbin_from_mean_and_std(0.6 * length as f64, 0.15 * length as f64)
                .unwrap();

        // Main generation step
        sample_multiple(amount, &rand::distr::Bernoulli::new(forward_seq_ratio).unwrap(), rng)
            .into_iter()
            .map(|is_forward| {
                let window_range: std::ops::Range<usize> = {
                    // generate window size
                    let window_size = w_size_distr.sample_single(rng).clamp(self.min_window_size, length);
                    // generate window starting index
                    let idx = rng.random_range(0..=length - window_size);
                    idx..idx + window_size
                };

                let mut seq = if is_forward {
                    core_seq[window_range.clone()].to_owned()
                } else {
                    core_seq_rc[window_range.clone()].to_owned()
                };

                // generate the block sizes of the quality scores
                let mut sum: usize = 0;
                let mut block_sizes: Vec<usize> = Vec::new();
                while sum < length {
                    let block_size = self.q_score_block_size_distr.sample_single(rng);
                    sum += block_size;
                    block_sizes.push(block_size);
                }
                block_sizes.last_mut().unwrap().sub_assign(sum - length);

                let mut q_scores: Vec<u8> = Vec::new();
                for block_size in block_sizes.into_iter() {
                    let central_q_score = self.q_score_distr.sample_single(rng);
                    q_scores.extend(
                        Self::Q_SCORE_INNER_BLOCK_OFFSET_DISTR
                            .sample_multiple(block_size, rng)
                            .into_iter()
                            .map(|val| {
                                ((val + central_q_score).saturating_sub(Self::ADJUSTMENT) as u8)
                                    .clamp(MIN_Q_SCORE, MAX_Q_SCORE)
                            }),
                    )
                }

                let mut nb_indel_errors = 0;
                let mut nb_insertion_errors = 0;
                let mut nb_snp_errors = 0;

                enum ErrorType {
                    Indel,
                    Insertion,
                    Snp,
                }

                let (indel_weight, insertion_weight, snp_weight) = self.indel_insertion_snp_error_weights;
                let error_type_to_weight = |etype: &ErrorType| -> f64 {
                    match etype {
                        ErrorType::Indel => indel_weight,
                        ErrorType::Insertion => insertion_weight,
                        ErrorType::Snp => snp_weight,
                    }
                };

                let mut idx = 0;
                while idx < seq.len() {
                    let p_error = myrio_core::constants::Q_TO_BP_CALL_ERROR_PROB_MAP[q_scores[idx] as usize];
                    if rng.sample(rand::distr::Bernoulli::new(p_error).unwrap()) {
                        // An error has "occurred", now we determine which type
                        match [ErrorType::Indel, ErrorType::Insertion, ErrorType::Snp]
                            .choose_weighted(rng, error_type_to_weight)
                            .unwrap()
                        {
                            ErrorType::Indel => {
                                seq.remove(idx..idx + 1);
                                q_scores.remove(idx);
                                nb_indel_errors += 1;
                            }
                            ErrorType::Insertion => {
                                seq.insert(
                                    idx,
                                    *[dna!["A"], dna!["C"], dna!["G"], dna!["T"]].choose(rng).unwrap(),
                                );
                                q_scores.insert(idx, q_scores[idx]);
                                nb_insertion_errors += 1;
                                idx += 1; // ignore the next one
                            }
                            ErrorType::Snp => {
                                let new_base = match seq.get(idx).unwrap() {
                                    Dna::A => [Dna::C, Dna::G, Dna::T],
                                    Dna::C => [Dna::A, Dna::G, Dna::T],
                                    Dna::G => [Dna::A, Dna::C, Dna::T],
                                    Dna::T => [Dna::A, Dna::C, Dna::G],
                                }
                                .choose(rng)
                                .unwrap()
                                .to_owned();

                                seq.set(idx, new_base);
                                q_scores[idx] = q_scores[idx].saturating_sub(2);
                                nb_snp_errors += 1;
                            }
                        }
                    }
                    idx += 1;
                }

                MyrSeq::new(
                    id.to_string(),
                    Some(format!("len={},wrange={window_range:?},ind={nb_indel_errors},ins={nb_insertion_errors},snp={nb_snp_errors}", seq.len())),
                    seq,
                    q_scores,
                )
            })
            .collect_vec()
    }
}

#[derive(Debug, Error)]
pub enum Error {
    #[error("Invalid coding to template strand ratio bounds, expected (a, b) such that 0 <= a <= b <= 1")]
    InvalidCTRBounds,
    #[error("Invalid error weights for ")]
    InvalidSequencingErrorWeights,
}
