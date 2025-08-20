// Modules
pub mod clustering;
pub mod constants;
pub mod data;
pub mod io;
pub mod simseq;
pub mod tax;

#[cfg(test)]
#[macro_export]
macro_rules! assert_float_eq {
    ($lhs: expr, $rhs: expr) => {
        let (a, b): (f64, f64) = ($lhs, $rhs);
        assert!((a - b).abs() < f64::EPSILON)
    };
}
