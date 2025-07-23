use std::collections::HashMap;

use argmin::{
    core::{CostFunction, Executor, Gradient},
    solver::{linesearch::MoreThuenteLineSearch, quasinewton::LBFGS},
};
use itertools::Itertools;
use myrio_core::MyrSeq;
use ndarray::{Array1, array};
use rand::SeedableRng;

use crate::simseq::generate_pseudo_amplicon;

const EPS_GRAD: f64 = 1e-2;

struct ClusterAlgorithm {
    k: usize,
    myrseqs: Vec<MyrSeq>,
}

impl ClusterAlgorithm {
    fn cluster(
        &self,
        t1_cutoff: f64,
        t2_cutoff: f64,
    ) -> Vec<Vec<MyrSeq>> {
        super::method_one(self.myrseqs.clone(), self.k, t1_cutoff, t2_cutoff, super::cosine_similarity)
            .unwrap()
    }
}

impl ClusterAlgorithm {
    fn cost_function(
        &self,
        t1_cutoff: f64,
        t2_cutoff: f64,
    ) -> f64 {
        let mut cost: f64 = 0.0;
        let clusters = self.cluster(t1_cutoff, t2_cutoff);

        for (idx, cluster) in clusters.into_iter().enumerate() {
            cost += idx as f64; // so that having many clusters get exceedingly worst
            let mut id_count_map: HashMap<&str, usize> = HashMap::new();
            for myrseq in cluster.iter() {
                let count_ref = id_count_map.entry(myrseq.id.as_str()).or_default();
                *count_ref += 1;
            }
            let counts = id_count_map.values().copied().collect_vec();
            let max_count = *counts.iter().max().unwrap();
            cost += ((counts.iter().sum::<usize>() - max_count) as f64) * 3.0;
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

pub fn run() -> anyhow::Result<()> {
    // Define cost function
    let cost = ClusterAlgorithm {
        myrseqs: MyrSeq::decode_vec_from_file("./ignore/argmin_myrseqs.bin").unwrap(),
        k: 5,
    };

    let init_param: Array1<f64> = array![0.66, 0.66];

    let init_grad = cost.gradient(&init_param)?;
    println!("Initial gradient: {init_grad:?}");

    let linesearch = MoreThuenteLineSearch::new().with_c(1e-4, 0.9)?;

    let solver = LBFGS::new(linesearch, 10); // 10 is memory size

    let res = Executor::new(cost, solver)
        .configure(|state| state.param(init_param.clone()).max_iters(100))
        .run()?;

    // Print result
    println!("{res}");
    Ok(())
}

pub fn _generate_and_save_random_sequences_for_cluster_testing() {
    let mut rng = rand::rngs::StdRng::from_os_rng();

    let myrseqs: Vec<MyrSeq> = [
        // ITS
        generate_pseudo_amplicon(&mut rng, 588, 150, "1"),
        // matK
        generate_pseudo_amplicon(&mut rng, 1500, 150, "2"),
        // rbcL
        generate_pseudo_amplicon(&mut rng, 1431, 150, "3"),
        // trnH-psbA
        generate_pseudo_amplicon(&mut rng, 1058, 150, "4"),
    ]
    .concat();

    MyrSeq::encode_vec_to_file("./ignore/argmin_myrseqs.bin", &myrseqs, 17).unwrap();
}
