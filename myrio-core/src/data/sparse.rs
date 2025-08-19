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

    pub fn merge_and_apply<F>(
        &self,
        other: &Self,
        op: F,
    ) -> Self
    where
        F: Fn((f64, f64)) -> f64,
    {
        let dim = self.dim.max(other.dim);
        let sval = op((self.sval, other.sval));
        let keys = self.merge_keys(other);
        let values = keys.iter().map(|&key| (self[key], other[key])).map(op).collect_vec();

        unsafe { Self::new_unchecked(keys, values, dim, sval) }
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
            let v = match a.get_unchecked(i).cmp(b.get_unchecked(j)) {
                Ordering::Less => {
                    i += 1;
                    *a.get_unchecked(i)
                }
                Ordering::Greater => {
                    j += 1;
                    *b.get_unchecked(j)
                }
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                    *a.get_unchecked(i)
                }
            };

            *ptr.add(written) = v;
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
