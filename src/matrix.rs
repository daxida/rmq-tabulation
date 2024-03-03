/// Matrix -- just indexing with i,j
pub struct Matrix {
    n: usize,
    table: Vec<usize>,
}

impl Matrix {
    pub fn new(n: usize) -> Matrix {
        let table = vec![0; n * n];
        Matrix { n, table }
    }

    fn flat_index(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < self.n && j < self.n);
        i * self.n + j
    }
}

impl std::ops::Index<(usize, usize)> for Matrix {
    type Output = usize;
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (i, j) = index;
        &self.table[self.flat_index(i, j)]
    }
}

impl std::ops::IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (i, j) = index;
        let flat_idx = self.flat_index(i, j);
        &mut self.table[flat_idx]
    }
}

/// Upper triangular table: Table for looking up at [i,j) (j > i) intervals.
#[derive(Clone)]
pub struct UTTable {
    n: usize,
    table: Vec<usize>,
}

impl UTTable {
    pub fn new(n: usize) -> UTTable {
        let table: Vec<usize> = vec![0; n * (n + 1) / 2];
        UTTable { n, table }
    }

    fn flat_index(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < self.n);
        debug_assert!(i < j && j <= self.n);
        let k = self.n - i - 1;
        k * (k + 1) / 2 + j - i - 1
    }
}

impl std::ops::Index<(usize, usize)> for UTTable {
    type Output = usize;
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (i, j) = index;
        &self.table[self.flat_index(i, j)]
    }
}

impl std::ops::IndexMut<(usize, usize)> for UTTable {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (i, j) = index;
        let flat_idx = self.flat_index(i, j);
        &mut self.table[flat_idx]
    }
}
