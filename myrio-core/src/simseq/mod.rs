// Modules
pub mod distr;

// Imports
use bio_seq::prelude::*;
use itertools::Itertools;
use rand::distr::Distribution;
use rand::seq::IndexedRandom;
use std::ops::SubAssign;

use crate::{
    MyrSeq,
    constants::{MAX_Q_SCORE, MIN_Q_SCORE},
};
use distr::{DiscreteDistribution, sample_multiple};

#[derive(Debug)]
pub struct Generator {
    pub coding_to_template_ratio_bounds: (f64, f64),
    pub q_score_distr: DiscreteDistribution,
    pub q_score_block_size_distr: DiscreteDistribution,
}

impl Default for Generator {
    fn default() -> Self {
        Self {
            coding_to_template_ratio_bounds: (0.35, 0.65),
            q_score_distr: DiscreteDistribution::new_cdf(&crate::constants::OBSERVED_Q_SCORE_CUMMUL_FREQ),
            q_score_block_size_distr: DiscreteDistribution::new_cdf(
                &crate::constants::OBSERVED_Q_SCORE_BLOCK_SIZE_CUMMUL_FREQ,
            ),
        }
    }
}

impl Generator {
    /// Represents a simple distribution of how q-scores are offset within a block of q-scores from the main q-score
    const Q_SCORE_INNER_BLOCK_OFFSET_DISTR: DiscreteDistribution =
        DiscreteDistribution::new_cdf(&[0.05, 0.15, 0.30, 0.70, 0.85, 0.95, 1.0]);

    /// `x = 0` is at index `3`
    const ADJUSTMENT: usize = 3;

    pub fn with_coding_to_template_ratio_bounds(
        self,
        bounds: (f64, f64),
    ) -> Self {
        Self { coding_to_template_ratio_bounds: bounds, ..self }
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

    pub fn generate_pseudo_amplicon(
        &self,
        length: usize,
        amount: usize,
        id: &str,
        rng: &mut impl rand::Rng,
    ) -> Vec<MyrSeq> {
        // We generate our core DNA forward amplicon sequence, and its reverse complement
        let core_seq: Seq<Dna> = {
            let mut buffer = vec![0_u8; ((usize::BITS as usize) >> 3) * length];
            rng.fill_bytes(&mut buffer);
            Seq::from_raw(length, bytemuck::cast_slice(&buffer)).unwrap()
        };
        let core_seq_rc = core_seq.to_revcomp();

        // We generate the ratio of sequences corresponding to the coding strand and the template (reverse complement) strand
        let forward_seq_ratio = rand::distr::uniform::Uniform::new(
            self.coding_to_template_ratio_bounds.0,
            self.coding_to_template_ratio_bounds.1,
        )
        .unwrap()
        .sample(&mut *rng);

        // Main generation step
        sample_multiple(amount, &rand::distr::Bernoulli::new(forward_seq_ratio).unwrap(), rng)
            .into_iter()
            .map(|is_forward| {
                let mut seq = if is_forward { core_seq.clone() } else { core_seq_rc.clone() };

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

                for (idx, q_score) in q_scores.iter().enumerate() {
                    let p_error = crate::constants::Q_TO_BP_CALL_ERROR_PROB_MAP[*q_score as usize];
                    if rng.sample(rand::distr::Bernoulli::new(p_error).unwrap()) {
                        seq.set(idx, *[Dna::A, Dna::C, Dna::G, Dna::T].choose(rng).unwrap())
                    }
                }

                MyrSeq::new(id.to_string(), None, seq, q_scores)
            })
            .collect_vec()
    }
}
