use itertools::Itertools;

#[derive(Debug, Clone, Default)]
pub struct SparseFloatVec {
    keys: Vec<usize>,
    values: Vec<f64>,
}

impl SparseFloatVec {
    pub fn new() -> Self {
        Self { keys: Vec::new(), values: Vec::new() }
    }

    pub fn get(
        &self,
        idx: usize,
    ) -> Option<&f64> {
        match self.keys.binary_search(&idx) {
            Ok(inner_idx) => Some(&self.values[inner_idx]),
            Err(_) => None,
        }
    }

    pub fn get_or_insert_mut(
        &mut self,
        idx: usize,
    ) -> &mut f64 {
        match self.keys.binary_search(&idx) {
            Ok(inner_idx) => &mut self.values[inner_idx],
            Err(inner_idx) => {
                self.keys.insert(inner_idx, idx);
                self.values.insert(inner_idx, f64::default());
                &mut self.values[inner_idx]
            }
        }
    }

    pub fn insert(
        &mut self,
        idx: usize,
        val: f64,
    ) -> Option<f64> {
        match self.keys.binary_search(&idx) {
            Ok(_) => Some(val),
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

    pub fn len(&self) -> usize {
        self.keys.len()
    }
}

impl core::ops::Index<usize> for SparseFloatVec {
    type Output = f64;

    fn index(
        &self,
        index: usize,
    ) -> &Self::Output {
        self.get(index).unwrap_or(&0.0)
    }
}
