// Imports
use core::{
    cmp::Ordering,
    ops::{AddAssign, DivAssign, MulAssign, SubAssign},
};
use std::{fmt::Debug, ops::Neg};

use bincode::{BorrowDecode, Decode, Encode};
use itertools::Itertools;
use myrio_proc::impl_ops_for_svec;

/// A sparse vector containing floats, similar to `HashMap<usize,Float>` but without the need for hashing
pub type Float = f64;
pub type SFVec = SparseVec<Float>;

#[derive(Debug, Clone, PartialEq)]
pub struct SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd,
{
    keys: Vec<usize>,
    values: Vec<T>,
    dim: usize,
    sval: T,
}

impl<T> SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd,
{
    pub fn new(dim: usize) -> Self {
        Self { keys: Vec::new(), values: Vec::new(), dim, sval: T::default() }
    }

    pub fn new_with_sparse_value(
        dim: usize,
        sval: T,
    ) -> Self {
        Self { keys: Vec::new(), values: Vec::new(), dim, sval }
    }

    /// # Safety
    /// `keys` and `values` must share the same length
    pub unsafe fn new_unchecked(
        keys: Vec<usize>,
        values: Vec<T>,
        dim: usize,
        sval: T,
    ) -> Self {
        Self { keys, values, dim, sval }
    }

    /// Returns the theoretical length of the array
    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn sval(&self) -> T {
        self.sval
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
    ) -> Option<&T> {
        match self.keys.binary_search(&idx) {
            Ok(inner_idx) => Some(unsafe { self.values.get_unchecked(inner_idx) }),
            Err(_) => None,
        }
    }

    pub fn get_or_insert_mut(
        &mut self,
        idx: usize,
    ) -> &mut T {
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
        val: T,
    ) -> Option<T> {
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

    pub fn iter(&self) -> impl Iterator<Item = (&usize, &T)> {
        self.keys.iter().zip_eq(self.values.iter())
    }

    pub fn keys(&self) -> impl Iterator<Item = &usize> {
        self.keys.iter()
    }

    pub fn values(&self) -> impl Iterator<Item = &T> {
        self.values.iter()
    }

    pub fn values_mut(&mut self) -> impl Iterator<Item = &mut T> {
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
        F: Fn(&mut T),
    {
        op(&mut self.sval);
        self.values.iter_mut().for_each(op);
    }

    pub fn merge_and_apply<F, O>(
        &self,
        other: &Self,
        op: F,
    ) -> SparseVec<O>
    where
        F: Fn(T, T) -> O,
        O: Debug + Clone + Copy + Default + PartialOrd,
    {
        let ak: &[usize] = &self.keys;
        let av: &[T] = &self.values;
        let asval: T = self.sval;
        let bk: &[usize] = &other.keys;
        let bv: &[T] = &other.values;
        let bsval: T = other.sval;

        let capacity: usize = ak.len() + bk.len();
        let mut ck: Vec<usize> = Vec::with_capacity(capacity);
        let mut cv: Vec<O> = Vec::with_capacity(capacity);

        unsafe {
            let ck_ptr: *mut usize = ck.as_mut_ptr();
            let cv_ptr: *mut O = cv.as_mut_ptr();
            let mut written = 0;

            let mut i = 0;
            let mut j = 0;

            while i < ak.len() && j < bk.len() {
                let ak_i = ak.get_unchecked(i);
                let av_i = av.get_unchecked(i);
                let bk_j = bk.get_unchecked(j);
                let bv_j = bv.get_unchecked(j);

                let (c_key, c_val): (usize, O) = match ak_i.cmp(bk_j) {
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

        unsafe { SparseVec::<O>::new_unchecked(ck, cv, dim, sval) }
    }
}

impl<T> SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd + AddAssign,
{
    pub fn from_unsorted_pairs(
        mut pairs: Vec<(usize, T)>,
        dim: usize,
        sval: T,
    ) -> Self {
        pairs.sort_unstable_by_key(|(key, _)| *key);
        let capacity = pairs.len();

        let mut keys: Vec<usize> = Vec::with_capacity(capacity);
        let mut values: Vec<T> = Vec::with_capacity(capacity);

        if capacity == 0 {
            return Self { keys, values, dim, sval };
        }

        unsafe {
            let k_ptr: *mut usize = keys.as_mut_ptr();
            let v_ptr: *mut T = values.as_mut_ptr();

            let mut idx: usize = 0;

            let (mut prev_key, mut val_to_write) = *pairs.get_unchecked(0);
            k_ptr.write(prev_key);

            pairs.into_iter().skip(1).for_each(|(key, val)| {
                if prev_key != key {
                    v_ptr.add(idx).write(val_to_write);
                    idx += 1;
                    k_ptr.add(idx).write(key);
                    val_to_write = val;
                    prev_key = key;
                } else {
                    val_to_write += val;
                }
            });
            v_ptr.add(idx).write(val_to_write);
            keys.set_len(idx + 1);
            //keys.shrink_to_fit();
            values.set_len(idx + 1);
            //values.shrink_to_fit();
        }

        Self { keys, values, dim, sval }
    }
}

impl SparseVec<Float> {
    pub fn sum(&self) -> Float {
        self.values.iter().sum1::<Float>().unwrap_or_default() + self.sval * self.nb_sparse() as Float
    }

    pub fn mean(&self) -> Float {
        self.sum() / self.dim as Float
    }

    pub fn norm_l2_squared(&self) -> Float {
        self.values().map(|v| v.powi(2)).sum1().unwrap_or(0.0) + self.sval.powi(2) * self.nb_sparse() as Float
    }

    pub fn norm_l2(&self) -> Float {
        self.norm_l2_squared().sqrt()
    }

    /// Alias for norm_l2, computes the magnitude (l2 norm) of the vector
    pub fn magnitude(&self) -> Float {
        self.norm_l2()
    }

    /// Normalize to a unit vector
    pub fn normalize_l2(&mut self) {
        let magn = self.magnitude();
        self.div_assign(magn);
    }

    pub fn normalize_l2_alt(mut self) -> Self {
        self.normalize_l2();
        self
    }

    /*
    pub fn clr_transform(&mut self) {
        self.apply_in_place(|x| x.add_assign(1.0));
        let ln_mean = self.mean().ln();
        self.apply_in_place(|x| *x = x.ln() - ln_mean);
    }

    /// https://journals.asm.org/doi/10.1128/msystems.00016-19#sec-4
    pub fn rclr_transform(&mut self) {
        let robust_ln_mean = self.values().product::<Float>().powf(1_Float / self.count() as Float).ln();
        assert!(robust_ln_mean.is_finite());
        self.values_mut().for_each(|x| *x = x.ln() - robust_ln_mean);
    }

    pub fn dist_l2(
        &self,
        other: &Self,
    ) -> Float {
        self.merge_and_apply(other, |x, y| (x - y).powi(2)).sum().sqrt()
    }
    */
}

impl<T> bincode::Encode for SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd + Encode,
{
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.keys, encoder)?;
        bincode::Encode::encode(&self.values, encoder)?;
        bincode::Encode::encode(&self.dim, encoder)?;
        bincode::Encode::encode(&self.sval, encoder)?;
        Ok(())
    }
}

impl<Context, T> bincode::Decode<Context> for SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd + Decode<Context> + 'static,
{
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(Self {
            keys: bincode::Decode::decode(decoder)?,
            values: bincode::Decode::decode(decoder)?,
            dim: bincode::Decode::decode(decoder)?,
            sval: bincode::Decode::decode(decoder)?,
        })
    }
}

impl<'de, Context, T> bincode::BorrowDecode<'de, Context> for SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd + BorrowDecode<'de, Context> + 'de,
{
    fn borrow_decode<D: bincode::de::BorrowDecoder<'de, Context = Context>>(
        decoder: &mut D
    ) -> Result<Self, bincode::error::DecodeError> {
        Ok(Self {
            keys: bincode::BorrowDecode::borrow_decode(decoder)?,
            values: bincode::BorrowDecode::borrow_decode(decoder)?,
            dim: bincode::BorrowDecode::borrow_decode(decoder)?,
            sval: bincode::BorrowDecode::borrow_decode(decoder)?,
        })
    }
}

impl<T> core::ops::Index<usize> for SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd,
{
    type Output = T;

    fn index(
        &self,
        index: usize,
    ) -> &Self::Output {
        self.get(index).unwrap_or(&self.sval)
    }
}

impl<T> core::ops::Neg for SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd + Neg<Output = T>,
{
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.values_mut().for_each(|v| *v = v.neg());
        self.sval = self.sval.neg();
        self
    }
}

impl<T> AsRef<SparseVec<T>> for SparseVec<T>
where
    T: Debug + Clone + Copy + Default + PartialOrd,
{
    fn as_ref(&self) -> &Self {
        self
    }
}

// run `cargo expand -p myrio-core data::sparse` to examine
impl_ops_for_svec!(Float Add);
impl_ops_for_svec!(Float Sub);
impl_ops_for_svec!(Float Mul);
impl_ops_for_svec!(Float Div);

impl<T> core::ops::Add<T> for SparseVec<Float>
where
    T: AsRef<SparseVec<Float>>,
{
    type Output = SparseVec<Float>;
    fn add(
        self,
        rhs: T,
    ) -> Self::Output {
        self.merge_and_apply(rhs.as_ref(), |x, y| x + y)
    }
}

impl<T> core::ops::Add<T> for &SparseVec<Float>
where
    T: AsRef<SparseVec<Float>>,
{
    type Output = SparseVec<Float>;
    fn add(
        self,
        rhs: T,
    ) -> Self::Output {
        self.merge_and_apply(rhs.as_ref(), |x, y| x + y)
    }
}

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

    use super::*;
    use crate::assert_float_eq;

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
        assert_eq!(sfvec.sval, Float::NEG_INFINITY)
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
        assert_float_eq!(sfvec.norm_l2(), 80_f64.sqrt())
    }
}
