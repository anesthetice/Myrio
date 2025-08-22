// Modules
pub mod clustering;
pub mod constants;
pub mod data;
pub mod io;
pub mod similarity;
pub mod simseq;
pub mod tax;

#[cfg(test)]
#[macro_export]
macro_rules! assert_float_eq {
    ($lhs: expr, $rhs: expr) => {
        assert!(($lhs - $rhs).abs() < f64::EPSILON)
    };
}
