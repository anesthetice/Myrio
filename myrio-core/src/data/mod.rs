// Modules
mod myrseq;
mod sparse;

// Re-exports
pub use myrseq::MyrSeq;
pub use sparse::{SFVec, SparseFloatVec};

pub type DFArray = ndarray::Array1<f64>;
