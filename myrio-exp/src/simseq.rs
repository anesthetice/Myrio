// Imports
use bio_seq::prelude::*;
use myrio_core::MyrSeq;
use rand_distr::{Distribution, Poisson};

pub fn generate_pseudo_amplicon(
    rng: &mut impl rand::Rng,
    length: usize,
    amount: usize,
) -> Vec<MyrSeq> {
    let core_seq: Seq<Dna> = {
        let mut buffer = vec![0_u8; ((usize::BITS as usize) >> 3) * length];
        rng.fill_bytes(&mut buffer);
        Seq::from_raw(length, bytemuck::cast_slice(&buffer)).unwrap()
    };

    let core_seq_rc = core_seq.to_revcomp();

    let Q_distrib = Poisson::new(2.0);
    let b = Q_distrib.unwrap().sample(rng) as usize;

    println!("{}", core_seq);

    Vec::new()
}
