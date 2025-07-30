use std::collections::HashMap;

use argmin::{
    core::{CostFunction, Executor, Gradient},
    solver::{linesearch::MoreThuenteLineSearch, quasinewton::LBFGS},
};
use itertools::Itertools;
use myrio_core::{MyrSeq, simseq::Generator};
use ndarray::{Array1, array};
use rand::SeedableRng;
use rayon::iter::ParallelIterator;
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator};

use crate::clustering::{ClusterMethod, SimFunc};

const EPS_GRAD: f64 = 2e-2;

struct ClusterAlgorithm {
    k: usize,
    myrseqs: Vec<MyrSeq>,
    method: ClusterMethod,
    sim_func: SimFunc,
}

impl ClusterAlgorithm {
    fn cluster(
        &self,
        t1_cutoff: f64,
        t2_cutoff: f64,
    ) -> Vec<Vec<MyrSeq>> {
        (self.method)(self.myrseqs.clone(), self.k, t1_cutoff, t2_cutoff, self.sim_func).unwrap()
    }
}

impl ClusterAlgorithm {
    fn cost_function(
        &self,
        t1_cutoff: f64,
        t2_cutoff: f64,
    ) -> f64 {
        let mut cost: f64 = 0.0;

        if !(0_f64..1_f64).contains(&t1_cutoff) || !(0_f64..1_f64).contains(&t2_cutoff) {
            return 1e10; // f64::MAX breaks `MoreThuenteLineSearch`
        }

        let clusters = self.cluster(t1_cutoff, t2_cutoff);

        for (idx, cluster) in clusters.into_iter().enumerate() {
            cost += idx as f64; // so that having many clusters gets exceedingly worst
            let mut id_count_map: HashMap<&str, usize> = HashMap::new();
            for myrseq in cluster.iter() {
                let count_ref = id_count_map.entry(myrseq.id.as_str()).or_default();
                *count_ref += 1;
            }
            let counts = id_count_map.values().copied().collect_vec();
            let max_count = *counts.iter().max().unwrap();
            let diff = (counts.iter().sum::<usize>() - max_count) as f64;
            cost += diff * diff * 2.0;
        }
        cost
    }
}

impl CostFunction for ClusterAlgorithm {
    type Output = f64;
    type Param = Array1<f64>;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, anyhow::Error> {
        let t1_cutoff = param[0];
        let t2_cutoff = param[1];
        Ok(self.cost_function(t1_cutoff, t2_cutoff))
    }
}

impl Gradient for ClusterAlgorithm {
    type Gradient = Array1<f64>;
    type Param = Array1<f64>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, anyhow::Error> {
        let t1 = param[0];
        let t2 = param[1];
        let base = self.cost_function(t1, t2);
        Ok(array![
            (self.cost_function(t1 + EPS_GRAD, t2) - base) / EPS_GRAD,
            (self.cost_function(t1, t2 + EPS_GRAD) - base) / EPS_GRAD,
        ])
    }
}

pub fn run_all(k: usize) -> anyhow::Result<()> {
    let myrseqs = MyrSeq::decode_vec_from_file("./ignore/argmin_myrseqs.bin")?;

    let costs = [
        ClusterAlgorithm {
            k,
            myrseqs: myrseqs.clone(),
            method: super::method_one,
            sim_func: super::cosine_similarity,
        },
        ClusterAlgorithm {
            k,
            myrseqs: myrseqs.clone(),
            method: super::method_one,
            sim_func: super::overlap_similarity,
        },
        ClusterAlgorithm {
            k,
            myrseqs: myrseqs.clone(),
            method: super::method_two,
            sim_func: super::cosine_similarity,
        },
        ClusterAlgorithm {
            k,
            myrseqs: myrseqs.clone(),
            method: super::method_two,
            sim_func: super::overlap_similarity,
        },
    ];

    let output = costs
        .into_par_iter()
        .map(|cost| {
            let init_param: Array1<f64> = array![0.7, 0.7];

            let init_grad = cost.gradient(&init_param)?;
            println!("Initial gradient: {init_grad:?}");

            let linesearch = MoreThuenteLineSearch::new().with_c(1e-4, 0.9)?;

            let solver = LBFGS::new(linesearch, 6);

            let res = Executor::new(cost, solver)
                .configure(|state| state.param(init_param.clone()).max_iters(100))
                .run()?;

            // Print result
            Ok(res.to_string())
        })
        .collect::<anyhow::Result<Vec<String>>>()?;

    println!("{}", output.iter().join("\n"));
    Ok(())
}

pub fn _generate_and_save_random_sequences_for_cluster_testing() {
    let mut rng = rand::rngs::StdRng::from_os_rng();
    let generator = Generator::default();

    let myrseqs: Vec<MyrSeq> = [
        // ITS
        generator.generate_pseudo_amplicon(588, 150, "1", &mut rng),
        // matK
        generator.generate_pseudo_amplicon(1500, 150, "2", &mut rng),
        // rbcL
        generator.generate_pseudo_amplicon(1431, 150, "3", &mut rng),
        // trnH-psbA
        generator.generate_pseudo_amplicon(1058, 150, "4", &mut rng),
    ]
    .concat();

    MyrSeq::encode_vec_to_file("./ignore/argmin_myrseqs.bin", &myrseqs, 17).unwrap();
}
