use itertools::Itertools;
// Imports
use myrio_core::{
    data::MyrSeq,
    simseq::{Generator, distr::DiscreteDistribution},
};
use rand::SeedableRng;

#[derive(Debug, Clone, Copy)]
#[repr(usize)]
enum Quality {
    Low = 8,
    Medium = 14,
    High = 24,
}

#[derive(Debug, Clone, Copy)]
#[repr(usize)]
enum Length {
    Short = 300,
    Medium = 900,
    Long = 2000,
}

const NB_MYRSEQ_PER_CLUSTER: usize = 100;
const QUALITIES: [Quality; 3] = [Quality::Low, Quality::Medium, Quality::High];
const LENGTHS: [Length; 3] = [Length::Short, Length::Medium, Length::Long];
const NB_FAMILIES: [usize; 3] = [2, 5, 10];

pub fn generate_testset(state_seed: u64) {
    fn generate_myrseqs(
        quality: Quality,
        length: Length,
        nb_families: usize,
        rng: &mut impl rand::Rng,
    ) -> Vec<MyrSeq> {
        let generator = Generator::default().with_q_score_distr(
            DiscreteDistribution::new_nbin_from_mean_and_std(
                quality as usize as f64,
                quality as usize as f64 * 0.4,
            )
            .unwrap(),
        );
        let mut myrseqs: Vec<MyrSeq> = Vec::new();
        for idx in 0..nb_families {
            myrseqs.extend(generator.generate_pseudo_amplicon(
                length as usize,
                NB_MYRSEQ_PER_CLUSTER,
                idx.to_string().as_str(),
                rng,
            ))
        }
        myrseqs
    }

    let mut rng = rand::rngs::StdRng::seed_from_u64(state_seed);

    for ((quality, length), nb_families) in
        QUALITIES.into_iter().cartesian_product(LENGTHS).cartesian_product(NB_FAMILIES)
    {
        let myrseqs = generate_myrseqs(quality, length, nb_families, &mut rng);
    }
}
