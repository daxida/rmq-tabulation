#[inline(always)]
fn flog2(v: usize) -> usize {
    v.ilog2() as usize
}

pub trait Rmq {
    fn rmq(&self, i: usize, j: usize) -> Option<usize>;
}

pub struct Sparse {
    lgn: usize,
    table: Vec<usize>,
}

impl Sparse {
    pub fn new(array: &[usize]) -> Self {
        let n = array.len();
        let lgn = flog2(n) + 1;
        let mut table: Vec<usize> = vec![0; lgn * n];

        for j in 0..n {
            table[j * lgn] = array[j];
        }

        for i in 1..lgn {
            for j in 0..=(n - (1 << i)) {
                table[j * lgn + i] = std::cmp::min(
                    table[j * lgn + i - 1],
                    table[(j + (1 << (i - 1))) * lgn + i - 1],
                );
            }
        }

        Sparse { lgn, table }
    }
}

impl Rmq for Sparse {
    fn rmq(&self, i: usize, j: usize) -> Option<usize> {
        if i >= j {
            return None;
        }

        let k = flog2(j - i);
        Some(std::cmp::min(
            self.table[i * self.lgn + k],
            self.table[(j - (1 << k)) * self.lgn + k],
        ))
    }
}

/// 2D Array class
struct Matrix {
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

/// Fully tabulating the answer to all queries with <O(n²),O(1)>
#[derive(Clone)]
struct TabulatedQuery {
    tbl: UTTable,
}

impl Rmq for TabulatedQuery {
    fn rmq(&self, i: usize, j: usize) -> Option<usize> {
        if i < j {
            Some(self.tbl[(i, j)])
        } else {
            None
        }
    }
}

impl TabulatedQuery {
    pub fn new(x: &[usize]) -> Self {
        let mut tbl = UTTable::new(x.len());
        for i in 0..x.len() {
            tbl[(i, i + 1)] = i;
        }
        for i in 0..x.len() - 1 {
            for j in i + 2..x.len() + 1 {
                // DP: min val in [i,j) is either min in [i,j-1) or [j-1,j)
                tbl[(i, j)] = if x[tbl[(i, j - 1)]] <= x[j - 1] {
                    tbl[(i, j - 1)]
                } else {
                    j - 1
                };
            }
        }
        TabulatedQuery { tbl }
    }
}

#[derive(Clone, Copy)]
struct BlockSize(pub usize);

#[derive(PartialEq, PartialOrd)]
struct BlockIdx(pub usize);

fn block_size(n: usize) -> BlockSize {
    BlockSize(flog2(n) + 1)
}

fn round_down(i: usize, bs: BlockSize) -> (BlockIdx, usize) {
    let BlockSize(bs) = bs;
    let r = i / bs;
    (BlockIdx(r), r * bs)
}
fn round_up(i: usize, bs: BlockSize) -> (BlockIdx, usize) {
    let BlockSize(bs) = bs;
    let r = (i + bs - 1) / bs;
    (BlockIdx(r), r * bs)
}

/// Reduce an array x to the smallest value in each block (of size block_size)
/// and the index in the original array that this minimal value sits at.
fn reduce_array(x: &[usize], block_size: BlockSize) -> Vec<usize> {
    let BlockSize(bs) = block_size;
    x.chunks(bs)
        .map(|block| *block.iter().min().unwrap())
        .collect()
}

/// Build a table of Ballot numbers B_pq from p=q=0 to p=q=b.
/// B_bb is the total number of binary trees with b leaves.
fn tabulate_ballot_numbers(b: usize) -> Matrix {
    let mut ballot = Matrix::new(b + 1);
    for q in 0..=b {
        ballot[(0, q)] = 1
    }
    for q in 1..=b {
        for p in 1..=q {
            ballot[(p, q)] = ballot[(p - 1, q)] + ballot[(p, q - 1)]
        }
    }
    ballot
}

/// Compute the block type number of a block.
/// The b argument is the true block size, but it can differ from block.len() for the last
/// block. When we process the last block, we fake push the missing elements, putting them
/// lower in the Cartesian tree than the real ones, so we still get the right RMQ.
fn block_type(block: &[usize], b: usize, stack: &mut [i64], ballot: &Matrix) -> usize {
    let mut num = 0;
    let mut top = 0;
    stack[top] = i64::MIN; // As close to -infinity as we get with this type...

    for (i, &v) in block.iter().enumerate() {
        let signed_v = v as i64;

        // Invariant: When i runs from zero to b, b-i is p in B[p,q].
        //            i-top is how the depth of the stack and b-(i-top)
        //            is then the q in B[p,q].
        let p = b - i;
        while stack[top] > signed_v {
            // Popping
            let q = b - (i - top);
            num += ballot[(p - 1, q)];
            top -= 1;
        }

        // Push...
        top += 1;
        stack[top] = signed_v;
    }

    num
}

/// Compute the block types and the tables for the blocks we observe in x.
fn tabulate_blocks(x: &[usize], b: usize) -> (Vec<usize>, Vec<Option<TabulatedQuery>>) {
    // We need to round up to get the number of blocks here.
    // The reduced array handles blocks 0, 1, ..., x.len()/b but we
    // might also have a block after it.
    let no_blocks = (x.len() + b - 1) / b;

    let ballot = tabulate_ballot_numbers(b);
    let mut stack: Vec<i64> = vec![i64::MIN; b + 1];

    let mut block_types = vec![0; no_blocks];
    let mut block_tables = vec![None; ballot[(b, b)]];
    for i in 0..no_blocks {
        let begin = i * b;
        // The last block might be smaller than a full block, but if we just
        // tabulate it anyway the missing values are virtually pushed and behave
        // like they are larger than the existing ones, giving us the right RMQ
        // results anyway (the true values are always smaller than the virtual ones).
        let end = std::cmp::min(x.len(), begin + b);
        let block = &x[begin..end];

        let bt = block_type(block, b, &mut stack, &ballot);
        block_types[i] = bt;
        if block_tables[bt].is_none() {
            block_tables[bt] = Some(TabulatedQuery::new(block));
        }
    }
    (block_types, block_tables)
}

pub struct Tabulation<'a> {
    x: &'a [usize],
    block_size: BlockSize,
    sparse: Sparse,
    block_types: Vec<usize>,
    block_tables: Vec<Option<TabulatedQuery>>,
}

impl<'a> Tabulation<'a> {
    pub fn new(x: &'a [usize]) -> Self {
        let n = x.len();
        let BlockSize(b) = block_size(n);

        // adjust block size; log(n) is too much so change it to a quarter.
        let b = std::cmp::max(4, b / 4); // I don't want too small blocks, so minimum is 4
        let block_size = BlockSize(b);

        let reduced_vals = reduce_array(x, block_size);
        let (block_types, block_tables) = tabulate_blocks(x, b);

        Tabulation {
            x,
            block_size,
            sparse: Sparse::new(&reduced_vals),
            block_types,
            block_tables,
        }
    }

    fn block_rmq(&self, i: usize, j: usize) -> Option<usize> {
        if j <= i {
            return None;
        }

        let BlockSize(bs) = self.block_size;
        let block_index = i / bs; // The index in the list of blocks
        let block_begin = block_index * bs; // The index the block starts at in x

        // Get the table for this block by looking up the block type...
        let btype = self.block_types[block_index];
        // ... and then the table from the block type.
        let tbl = self.block_tables[btype].as_ref().unwrap();

        // Get RMQ and adjust the index back up, so it is relative to the start of the block.
        let rmq_idx = tbl.rmq(i - block_begin, j - block_begin)? + block_begin;
        Some(self.x[rmq_idx])
    }
}

fn lift_op<T: Copy>(f: impl Fn(T, T) -> T) -> impl Fn(Option<T>, Option<T>) -> Option<T> {
    move |a, b| match (a, b) {
        (None, None) => None,
        (Some(_), None) => a,
        (None, Some(_)) => b,
        (Some(a), Some(b)) => Some(f(a, b)),
    }
}

impl<'a> Rmq for Tabulation<'a> {
    fn rmq(&self, i: usize, j: usize) -> Option<usize> {
        let BlockSize(bs) = self.block_size;
        let bi = BlockIdx(i / bs);
        // The block indices are not the same for the small tables and the sparse table.
        // For the sparse table we have to round up for i ...
        let (sparse_bi, ii) = round_up(i, BlockSize(bs));
        // ... but to get the block i is in, we need to round down.
        let (bj, jj) = round_down(j, BlockSize(bs));

        if bi < bj {
            let p1 = self.block_rmq(i, ii);
            let (BlockIdx(sparse_bi), BlockIdx(bj)) = (sparse_bi, bj);
            let p2 = self.sparse.rmq(sparse_bi, bj);
            let p3 = self.block_rmq(jj, j);
            let min = lift_op(std::cmp::min);
            Some(min(min(p1, p2), p3)?)
        } else {
            Some(self.block_rmq(i, j)?)
        }
    }
}

#[cfg(test)]
mod tests;
