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
use distr::{FloatDistribution, UsizeDistribution, sample_multiple};

#[derive(Debug)]
pub struct Generator {
    pub coding_to_template_ratio_bounds: (f64, f64),
    pub q_score_distr: UsizeDistribution,
    pub q_score_block_size_distr: UsizeDistribution,
    pub q_score_block_inner_std_dev: f64,
}

impl Default for Generator {
    fn default() -> Self {
        Self {
            coding_to_template_ratio_bounds: (0.3, 0.7),
            q_score_distr: UsizeDistribution::new_cdf(&crate::constants::Q_SCORE_CUMMUL_FREQ),
            q_score_block_size_distr: UsizeDistribution::new_cdf(
                &crate::constants::Q_SCORE_BLOCK_SIZE_CUMMUL_FREQ,
            ),
            q_score_block_inner_std_dev: 2.0,
        }
    }
}

impl Generator {
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
                        sample_multiple::<f64>(
                            block_size,
                            &rand_distr::Normal::new(
                                central_q_score as f64,
                                self.q_score_block_inner_std_dev,
                            )
                            .unwrap(),
                            rng,
                        )
                        .into_iter()
                        .map(|val| (val.round() as u8).clamp(MIN_Q_SCORE, MAX_Q_SCORE)),
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
