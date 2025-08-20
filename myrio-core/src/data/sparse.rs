// Imports
use core::{
    cmp::Ordering,
    ops::{AddAssign, DivAssign, MulAssign, SubAssign},
};

use itertools::Itertools;
use myrio_proc::impl_f64_ops_for_sfvec;

pub type SFVec = SparseFloatVec;

/// A sparse vector containing floats, similar to `HashMap<usize,f64>` but without the need for hashing
#[derive(Debug, Clone)]
pub struct SparseFloatVec {
    keys: Vec<usize>,
    values: Vec<f64>,
    dim: usize,
    /// Baseline sparse value, note that this does not act as an offset for 'real' values
    pub sval: f64,
}

impl SparseFloatVec {
    pub fn new(dim: usize) -> Self {
        Self { keys: Vec::new(), values: Vec::new(), dim, sval: 0.0 }
    }

    pub fn new_with_sparse_value(
        dim: usize,
        sval: f64,
    ) -> Self {
        Self { keys: Vec::new(), values: Vec::new(), dim, sval }
    }

    unsafe fn new_unchecked(
        keys: Vec<usize>,
        values: Vec<f64>,
        dim: usize,
        sval: f64,
    ) -> Self {
        Self { keys, values, dim, sval }
    }

    /// Returns the theoretical length of the array
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Returns the number of non-sparse values held
    pub fn count(&self) -> usize {
        self.keys.len()
    }

    pub fn nb_sparse(&self) -> usize {
        self.dim - self.count()
    }

    pub fn get(
        &self,
        idx: usize,
    ) -> Option<&f64> {
        match self.keys.binary_search(&idx) {
            Ok(inner_idx) => Some(unsafe { self.values.get_unchecked(inner_idx) }),
            Err(_) => None,
        }
    }

    pub fn get_or_insert_mut(
        &mut self,
        idx: usize,
    ) -> &mut f64 {
        if idx >= self.dim() {
            panic!("Specified index of {idx} exceeds SFVec dimension of {}", self.dim)
        }
        match self.keys.binary_search(&idx) {
            Ok(inner_idx) => unsafe { self.values.get_unchecked_mut(inner_idx) },
            Err(inner_idx) => {
                self.keys.insert(inner_idx, idx);
                self.values.insert(inner_idx, self.sval);
                unsafe { self.values.get_unchecked_mut(inner_idx) }
            }
        }
    }

    pub fn insert(
        &mut self,
        idx: usize,
        val: f64,
    ) -> Option<f64> {
        if idx >= self.dim() {
            panic!("Specified index of {idx} exceeds SFVec dimension of {}", self.dim)
        }
        match self.keys.binary_search(&idx) {
            Ok(inner_idx) => {
                self.values.push(val);
                Some(self.values.swap_remove(inner_idx))
            }
            Err(inner_idx) => {
                self.keys.insert(inner_idx, idx);
                self.values.insert(inner_idx, val);
                None
            }
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (&usize, &f64)> {
        self.keys.iter().zip_eq(self.values.iter())
    }

    pub fn keys(&self) -> impl Iterator<Item = &usize> {
        self.keys.iter()
    }

    pub fn values(&self) -> impl Iterator<Item = &f64> {
        self.values.iter()
    }

    pub fn values_mut(&mut self) -> impl Iterator<Item = &mut f64> {
        self.values.iter_mut()
    }

    pub fn merge_keys(
        &self,
        other: &Self,
    ) -> Vec<usize> {
        unsafe { merge_sorted_unique(&self.keys, &other.keys) }
    }

    pub fn apply_in_place<F>(
        &mut self,
        op: F,
    ) where
        F: Fn(&mut f64),
    {
        op(&mut self.sval);
        self.values.iter_mut().for_each(op);
    }

    pub fn merge_and_apply<F>(
        &self,
        other: &Self,
        op: F,
    ) -> Self
    where
        F: Fn(f64, f64) -> f64,
    {
        let ak: &[usize] = &self.keys;
        let av: &[f64] = &self.values;
        let asval: f64 = self.sval;
        let bk: &[usize] = &other.keys;
        let bv: &[f64] = &other.values;
        let bsval: f64 = other.sval;

        let capacity: usize = ak.len() + bk.len();
        let mut ck: Vec<usize> = Vec::with_capacity(capacity);
        let mut cv: Vec<f64> = Vec::with_capacity(capacity);

        unsafe {
            let ck_ptr: *mut usize = ck.as_mut_ptr();
            let cv_ptr: *mut f64 = cv.as_mut_ptr();
            let mut written = 0;

            let mut i = 0;
            let mut j = 0;

            while i < ak.len() && j < bk.len() {
                let ak_i = ak.get_unchecked(i);
                let av_i = av.get_unchecked(i);
                let bk_j = bk.get_unchecked(j);
                let bv_j = bv.get_unchecked(j);

                let (c_key, c_val): (usize, f64) = match ak_i.cmp(bk_j) {
                    Ordering::Less => {
                        i += 1;
                        (*ak_i, op(*av_i, bsval))
                    }
                    Ordering::Greater => {
                        j += 1;
                        (*bk_j, op(*bv_j, asval))
                    }
                    Ordering::Equal => {
                        i += 1;
                        j += 1;
                        (*ak_i, op(*av_i, *bv_j))
                    }
                };

                *ck_ptr.add(written) = c_key;
                *cv_ptr.add(written) = c_val;
                written += 1;
            }

            // drain remaining elements either in `a` or `b` (will never be both)
            while i < ak.len() {
                *ck_ptr.add(written) = *ak.get_unchecked(i);
                *cv_ptr.add(written) = op(*av.get_unchecked(i), bsval);
                written += 1;
                i += 1;
            }
            while j < bk.len() {
                *ck_ptr.add(written) = *bk.get_unchecked(j);
                *cv_ptr.add(written) = op(*bv.get_unchecked(j), asval);
                written += 1;
                j += 1;
            }

            // tell Vec how many elements we actually wrote
            ck.set_len(written);
            cv.set_len(written);
        }

        let dim = self.dim.max(other.dim);
        let sval = op(asval, bsval);

        unsafe { Self::new_unchecked(ck, cv, dim, sval) }
    }

    pub fn sum(&self) -> f64 {
        self.values().sum1().unwrap_or(0.0) + self.sval * self.nb_sparse() as f64
    }

    pub fn mean(&self) -> f64 {
        self.sum() / self.dim as f64
    }

    pub fn clr_transform(&mut self) {
        self.apply_in_place(|x| x.add_assign(1.0));
        let ln_mean = self.mean().ln();
        self.apply_in_place(|x| *x = x.ln() - ln_mean);
    }

    pub fn dist_l2(
        &self,
        other: &Self,
    ) -> f64 {
        self.merge_and_apply(other, |x, y| (x - y).powi(2)).sum().sqrt()
    }
}

impl core::ops::Index<usize> for SparseFloatVec {
    type Output = f64;

    fn index(
        &self,
        index: usize,
    ) -> &Self::Output {
        self.get(index).unwrap_or(&self.sval)
    }
}

impl core::ops::Neg for SparseFloatVec {
    type Output = SparseFloatVec;

    fn neg(mut self) -> Self::Output {
        self.values_mut().for_each(|v| *v = v.neg());
        self.sval = self.sval.neg();
        self
    }
}

// run `cargo expand -p myrio-core data::sparse` to examine
impl_f64_ops_for_sfvec!(Add);
impl_f64_ops_for_sfvec!(Sub);
impl_f64_ops_for_sfvec!(Mul);
impl_f64_ops_for_sfvec!(Div);

/// This function is only safe if both slices are sorted with unique elements
unsafe fn merge_sorted_unique(
    a: &[usize],
    b: &[usize],
) -> Vec<usize> {
    let mut out = Vec::with_capacity(a.len() + b.len());

    unsafe {
        let ptr: *mut usize = out.as_mut_ptr(); // raw pointer to output buffer
        let mut written = 0;

        let mut i = 0;
        let mut j = 0;

        while i < a.len() && j < b.len() {
            let a_i = a.get_unchecked(i);
            let b_j = b.get_unchecked(j);
            let v = match a_i.cmp(b_j) {
                Ordering::Less => {
                    i += 1;
                    a_i
                }
                Ordering::Greater => {
                    j += 1;
                    b_j
                }
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                    a_i
                }
            };

            *ptr.add(written) = *v;
            written += 1;
        }

        // drain remaining elements either in `a` or `b` (will never be both)
        while i < a.len() {
            *ptr.add(written) = *a.get_unchecked(i);
            written += 1;
            i += 1;
        }
        while j < b.len() {
            *ptr.add(written) = *b.get_unchecked(j);
            written += 1;
            j += 1;
        }

        // tell Vec how many elements we actually wrote
        out.set_len(written);
    }

    out
}

#[cfg(test)]
mod test {
    use std::f64::consts::E;

    use crate::assert_float_eq;

    use super::*;

    #[test]
    fn merge_sorted_unique_fn_test() {
        let a: [usize; 4] = [1, 2, 6, 7];
        let b: [usize; 8] = [2, 3, 4, 9, 10, 11, 12, 13];

        assert_eq!(unsafe { merge_sorted_unique(&a, &b) }, vec![1, 2, 3, 4, 6, 7, 9, 10, 11, 12, 13])
    }

    #[test]
    fn apply_in_place_method_test() {
        let mut sfvec = unsafe { SFVec::new_unchecked(vec![1, 2], vec![1.0, E], 2, 0.0) };
        sfvec.apply_in_place(|v| *v = v.ln());
        assert_float_eq!(sfvec[1], 0.0);
        assert_float_eq!(sfvec[2], 1.0);
        assert_eq!(sfvec.sval, f64::NEG_INFINITY)
    }

    #[test]
    fn merge_and_apply_method_test() {
        let a = unsafe { SFVec::new_unchecked(vec![1, 2, 3], vec![1.0, 2.0, 3.0], 3, 1.0) };
        let b = unsafe { SFVec::new_unchecked(vec![2, 3, 4], vec![2.0, 3.0, 4.0], 3, 2.0) };

        let c = a.merge_and_apply(&b, |a, b| a + b);

        assert_float_eq!(c[1], 3.0);
        assert_float_eq!(c[2], 4.0);
        assert_float_eq!(c[3], 6.0);
        assert_float_eq!(c[4], 5.0);
    }

    #[test]
    fn misc_method_tests() {
        let sfvec = unsafe { SFVec::new_unchecked(vec![1, 2], vec![1.0, 2.0], 5, 5.0) };

        assert_float_eq!(sfvec.sum(), 18.0);
        assert_float_eq!(sfvec.mean(), 18.0 / 5.0);
    }
}
