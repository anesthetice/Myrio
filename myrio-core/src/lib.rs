// Modules
pub mod clustering;
pub mod constants;
pub mod data;
pub mod io;
pub mod prelude;
pub mod similarity;
pub mod tax;
pub mod utils;

#[cfg(test)]
#[macro_export]
macro_rules! assert_float_eq {
    ($lhs: expr, $rhs: expr) => {
        assert!(($lhs - $rhs).abs() < f32::EPSILON)
    };
}
