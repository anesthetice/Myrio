use argmin::core::CostFunction;

enum SimilarityFunction {
    Cosine,
    Overlap,
}

struct ClusterAlgorithm {}

impl CostFunction for ClusterAlgorithm {
    type Param = (f64, f64, SimilarityFunction);
    type Output = f64;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, anyhow::Error> {
        unimplemented!()
    }
}

fn generate_random_seqs() {}
