use std::ops::SubAssign;

// Imports
use bio_seq::prelude::*;
use itertools::Itertools;
use myrio_core::MyrSeq;
use rand::seq::IndexedRandom;
use rand_distr::{Distribution, Poisson};

pub fn generate_pseudo_amplicon(
    mut rng: &mut impl rand::Rng,
    length: usize,
    amount: usize,
    id: &str,
) -> Vec<MyrSeq> {
    // We generate our core DNA forward amplicon sequence, and its reverse complement
    let core_seq: Seq<Dna> = {
        let mut buffer = vec![0_u8; ((usize::BITS as usize) >> 3) * length];
        rng.fill_bytes(&mut buffer);
        Seq::from_raw(length, bytemuck::cast_slice(&buffer)).unwrap()
    };
    let core_seq_rc = core_seq.to_revcomp();

    // We generate the ratio of forward to revcomp sequence
    let forward_seq_ratio = rand::distr::uniform::Uniform::new(0.3, 0.7).unwrap().sample(&mut *rng);

    // Main generation step
    gen_n_samples(rng, &rand::distr::Bernoulli::new(forward_seq_ratio).unwrap(), amount)
        .into_iter()
        .map(|is_forward| {
            let mut seq = if is_forward { core_seq.clone() } else { core_seq_rc.clone() };

            // generate the block sizes of the quality scores
            let mut sum: usize = 0;
            let mut block_sizes: Vec<usize> = Vec::new();
            while sum < length {
                let rand_unif_val: f64 = rng.sample(rand::distr::StandardUniform);
                let block_size = Q_SCORE_BLOCK_SIZE_CUMMUL_FREQ
                    .into_iter()
                    .find_position(|freq| *freq >= rand_unif_val)
                    .unwrap()
                    .0;
                block_sizes.push(block_size);
                sum += block_size;
            }
            block_sizes.last_mut().unwrap().sub_assign(sum - length);

            let mut q_scores: Vec<u8> = Vec::new();
            for block_size in block_sizes.into_iter() {
                let rand_unif_val: f64 = rng.sample(rand::distr::StandardUniform);
                let central_q_score =
                    Q_SCORE_CUMMUL_FREQ.into_iter().find_position(|freq| *freq >= rand_unif_val).unwrap().0
                        as f64;
                q_scores.extend(
                    gen_n_samples::<f64>(
                        rng,
                        &rand_distr::Normal::new(central_q_score, 2.0).unwrap(),
                        block_size,
                    )
                    .into_iter()
                    .map(|val| (val.round() as u8).clamp(MIN_Q_SCORE, MAX_Q_SCORE)),
                )
            }

            for (idx, q_score) in q_scores.iter().enumerate() {
                let p_error = myrio_core::constants::Q_TO_BP_CALL_ERROR_PROB_MAP[*q_score as usize];
                if rng.sample(rand::distr::Bernoulli::new(p_error).unwrap()) {
                    seq.set(idx, *[Dna::A, Dna::C, Dna::G, Dna::T].choose(rng).unwrap())
                }
            }

            MyrSeq::new(seq, q_scores, Some(id.to_string()), None)
        })
        .collect_vec()
}

fn gen_n_samples<T>(
    rng: &mut impl rand::Rng,
    distrib: &impl rand::distr::Distribution<T>,
    n: usize,
) -> Vec<T> {
    let mut output: Vec<T> = Vec::new();
    for _ in 0..n {
        output.push(rng.sample(distrib));
    }
    output
}

#[rustfmt::skip]
const Q_SCORE_CUMMUL_FREQ: [f64; 41] = [
    0.0, 2.68953393e-04, 5.69554715e-03, 2.43616857e-02,
    5.69442911e-02, 9.83357123e-02, 1.43028724e-01, 1.85568252e-01,
    2.24820989e-01, 2.61233883e-01, 2.95625883e-01, 3.28616438e-01,
    3.60643874e-01, 3.92169862e-01, 4.23491155e-01, 4.54940339e-01,
    4.86660050e-01, 5.18859644e-01, 5.51574575e-01, 5.84902833e-01,
    6.18956183e-01, 6.53794544e-01, 6.89246839e-01, 7.25249907e-01,
    7.61609172e-01, 7.97842574e-01, 8.33153333e-01, 8.66356580e-01,
    8.96342126e-01, 9.22066069e-01, 9.43142934e-01, 9.59716815e-01,
    9.72279533e-01, 9.81483394e-01, 9.87986881e-01, 9.92487061e-01,
    9.95518922e-01, 9.97476939e-01, 9.98729194e-01, 9.99504143e-01,
    1.0,
];

#[rustfmt::skip]
const Q_SCORE_BLOCK_SIZE_CUMMUL_FREQ: [f64; 15] = [
    0.0, 0.34458653, 0.57879562, 0.73259496,
    0.83996226, 0.89905423, 0.93538613, 0.95836335,
    0.97319795, 0.98278344, 0.98912095, 0.99347341,
    0.99649658, 0.99855931, 1.0
];

const MAX_Q_SCORE: u8 = 40;
const MIN_Q_SCORE: u8 = 1; // Technically 0, but we'll use 1 in this case as we never observe 0
