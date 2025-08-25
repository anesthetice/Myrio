// Imports

use itertools::Itertools;
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

impl From<&str> for Quality {
    fn from(value: &str) -> Self {
        match value {
            "Low" => Self::Low,
            "Medium" => Self::Medium,
            "High" => Self::High,
            _ => unreachable!(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
#[repr(usize)]
enum Length {
    Short = 300,
    Medium = 900,
    Long = 2000,
}

impl From<&str> for Length {
    fn from(value: &str) -> Self {
        match value {
            "Short" => Self::Short,
            "Medium" => Self::Medium,
            "Long" => Self::Long,
            _ => unreachable!(),
        }
    }
}

const NB_MYRSEQ_PER_CLUSTER: usize = 100;
const QUALITIES: [Quality; 3] = [Quality::Low, Quality::Medium, Quality::High];
const LENGTHS: [Length; 3] = [Length::Short, Length::Medium, Length::Long];
const NB_FAMILIES: [usize; 3] = [2, 5, 10];

pub fn generate_testset(state_seed: u64) -> anyhow::Result<()> {
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
        let filepath = std::path::PathBuf::from(format!(
            "./ignore/testset/len={length:?}_qual={quality:?}_fam={nb_families}"
        ));
        println!("\"{}\",", filepath.display());
        MyrSeq::encode_vec_to_file(filepath, &myrseqs, 20)?;
    }
    Ok(())
}

pub fn load_testset() -> anyhow::Result<Vec<(Vec<MyrSeq>, usize, String)>> {
    const FILEPATHS: [&str; 27] = [
        "./ignore/testset/len=Short_qual=Low_fam=2",
        "./ignore/testset/len=Short_qual=Low_fam=5",
        "./ignore/testset/len=Short_qual=Low_fam=10",
        "./ignore/testset/len=Medium_qual=Low_fam=2",
        "./ignore/testset/len=Medium_qual=Low_fam=5",
        "./ignore/testset/len=Medium_qual=Low_fam=10",
        "./ignore/testset/len=Long_qual=Low_fam=2",
        "./ignore/testset/len=Long_qual=Low_fam=5",
        "./ignore/testset/len=Long_qual=Low_fam=10",
        "./ignore/testset/len=Short_qual=Medium_fam=2",
        "./ignore/testset/len=Short_qual=Medium_fam=5",
        "./ignore/testset/len=Short_qual=Medium_fam=10",
        "./ignore/testset/len=Medium_qual=Medium_fam=2",
        "./ignore/testset/len=Medium_qual=Medium_fam=5",
        "./ignore/testset/len=Medium_qual=Medium_fam=10",
        "./ignore/testset/len=Long_qual=Medium_fam=2",
        "./ignore/testset/len=Long_qual=Medium_fam=5",
        "./ignore/testset/len=Long_qual=Medium_fam=10",
        "./ignore/testset/len=Short_qual=High_fam=2",
        "./ignore/testset/len=Short_qual=High_fam=5",
        "./ignore/testset/len=Short_qual=High_fam=10",
        "./ignore/testset/len=Medium_qual=High_fam=2",
        "./ignore/testset/len=Medium_qual=High_fam=5",
        "./ignore/testset/len=Medium_qual=High_fam=10",
        "./ignore/testset/len=Long_qual=High_fam=2",
        "./ignore/testset/len=Long_qual=High_fam=5",
        "./ignore/testset/len=Long_qual=High_fam=10",
    ];

    FILEPATHS
        .into_iter()
        .map(|fp| {
            let myrseqs = MyrSeq::decode_vec_from_file(fp)?;
            let length = Length::from(fp.split_once("len=").unwrap().1.split_once("_").unwrap().0);
            let quality = Quality::from(fp.split_once("qual=").unwrap().1.split_once("_").unwrap().0);
            let nb_families = fp.split_once("fam=").unwrap().1.parse()?;
            Ok((myrseqs, nb_families, format!("length: {length:?}, quality: {quality:?}, nf: {nb_families}")))
        })
        .collect()
}
