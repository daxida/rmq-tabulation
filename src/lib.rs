#[inline(always)]
fn flog2(v: usize) -> usize {
    v.ilog2() as usize
}

pub trait RMQ {
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

impl RMQ for Sparse {
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

mod matrix;
use matrix::Matrix;

/// An index i together with its value x[i].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Point(usize, usize);
impl Point {
    pub fn new(i: usize, x: &[usize]) -> Point {
        Point(i, x[i])
    }
    /// Get Some(Point(i,x[i])) if the index is valid or None if not
    pub fn get(i: Option<usize>, x: &[usize]) -> Option<Point> {
        Some(Point(i?, *x.get(i?)?))
    }
}

impl std::cmp::PartialOrd for Point {
    fn partial_cmp(&self, other: &Point) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl std::cmp::Ord for Point {
    fn cmp(&self, other: &Point) -> std::cmp::Ordering {
        self.1.cmp(&other.1).then_with(|| self.0.cmp(&other.0))
    }
}

use matrix::UTTable;

/// Fully tabulating the answer to all queries with
/// <O(n²),O(1)> running times
#[derive(Clone)]
pub struct TabulatedQuery {
    tbl: UTTable,
}

impl RMQ for TabulatedQuery {
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
                // Dynamic programming:
                // Min val in [i,j) is either min in [i,j-1) or [j-1,j)
                let left = Point::new(tbl[(i, j - 1)], x);
                let current = Point::new(j - 1, x);
                tbl[(i, j)] = std::cmp::min(left, current).0
            }
        }
        TabulatedQuery { tbl }
    }
}

/// Storing values for all (i,i+2^k) indices and only those
pub mod powers {
    /// Type for powers of two, 2^k. Contains k, but wrapped in
    /// a type so we don't confuse log-space with linear space.
    #[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
    pub struct Pow(pub usize);

    impl Pow {
        /// for a power Pow(k) get 2^k.
        #[inline]
        pub fn value(&self) -> usize {
            1 << self.0
        }
    }

    /// Get k such that 2**k is j rounded down to the
    /// nearest power of 2.
    /// j=1=2^0 => 0
    /// j=2=2^1 => 1
    /// j=3=2^1+1 => 1
    /// j=4=2^2 => 2
    /// and so on.
    pub fn log2_down(j: usize) -> Pow {
        assert!(j != 0); // not defined for zero

        // Rounded down means finding the index of the first
        // 1 in the bit-pattern. If j = 00010101110
        // then 00010000000 (only first bit) is the closest
        // power of two, and we want the position of that bit.
        // j.leading_zeros() counts the number of leading zeros
        // and we get the index by subtracting this
        // from the total number of bits minus one.
        Pow((usize::BITS - j.leading_zeros() - 1) as usize)
        // usize::BITS and j.leading_zeros() will be u32, so
        // we cast the result back to usize.
    }

    pub fn power_of_two(x: usize) -> bool {
        (x == 0) || ((x & (x - 1)) == 0)
    }

    /// For n, get (rounded up) log2(n).
    pub fn log2_up(n: usize) -> Pow {
        // log_table_size(n) with n=2^k+m will always give us 2^{k+1},
        // whether m is zero or not. We want 2^{k+1} when m > 0 and 2^k
        // when m is zero, i.e. when n is a power of two.
        // So we should subtract one from the exponent if n is a power of two.
        let Pow(k) = log_table_size(n);
        Pow(k - power_of_two(n) as usize)
    }

    /// We always have to add one to the exponent, because in log-space
    /// we are working with 1-indexed (0-indexed in log-space) values,
    /// so to have a table that can handle maximum value k, we need k+1
    /// entires. That is what this function gives us.
    pub fn log_table_size(n: usize) -> Pow {
        let Pow(k) = log2_down(n);
        Pow(k + 1)
    }

    /// From range [i,j), get values (k,j-2^k) where k is the offset
    /// into the TwoD table to look up the value for [i,i+2^k) and [j-2^k,j)
    /// from which we can get the RMQ.
    pub fn power_index(i: usize, j: usize) -> ((usize, Pow), (usize, Pow)) {
        let powk = log2_down(j - i);
        ((i, powk), (j - powk.value(), powk))
    }

    /// A rather simple 2D array made from vectors of vectors.
    /// There are better solutions, but I can implement those later
    /// with the same interface.
    pub struct Powers {
        table: Vec<Vec<usize>>,
    }

    impl Powers {
        pub fn new(n: usize) -> Powers {
            let Pow(logn) = log_table_size(n);
            let table = vec![vec![0; logn]; n];
            Powers { table }
        }
    }

    impl std::ops::Index<(usize, Pow)> for Powers {
        type Output = usize;
        fn index(&self, index: (usize, Pow)) -> &Self::Output {
            let (i, Pow(k)) = index;
            &self.table[i][k]
        }
    }

    impl std::ops::IndexMut<(usize, Pow)> for Powers {
        fn index_mut(&mut self, index: (usize, Pow)) -> &mut Self::Output {
            let (i, Pow(k)) = index;
            &mut self.table[i][k]
        }
    }
}

#[derive(Clone, Copy)]
pub struct BlockSize(pub usize);

#[derive(PartialEq, PartialOrd)]
pub struct BlockIdx(pub usize);

pub fn block_size(n: usize) -> BlockSize {
    let powers::Pow(block_size) = powers::log2_up(n);
    BlockSize(block_size)
}

pub fn round_down(i: usize, bs: BlockSize) -> (BlockIdx, usize) {
    let BlockSize(bs) = bs;
    let r = i / bs;
    (BlockIdx(r), r * bs)
}
pub fn round_up(i: usize, bs: BlockSize) -> (BlockIdx, usize) {
    let BlockSize(bs) = bs;
    let r = (i + bs - 1) / bs;
    (BlockIdx(r), r * bs)
}

/// Reduce an array x to the smallest value in each block (of size block_size)
/// and the index in the original array that this minimal value sits at.
pub fn reduce_array(x: &[usize], block_size: BlockSize) -> (Vec<usize>, Vec<usize>) {
    let BlockSize(bs) = block_size;
    let mut indices: Vec<usize> = Vec::new();
    let mut values: Vec<usize> = Vec::new();
    let no_blocks = x.len() / bs;
    for block in 0..no_blocks {
        let begin = block * bs;
        let end = begin + bs;
        // naive rmq
        let y = &x[begin..end];
        let min_val = y.iter().min().unwrap();
        let pos = begin + y.iter().position(|a| a == min_val).unwrap();
        let Point(pos, val) = Point::new(pos, x);
        indices.push(pos);
        values.push(val);
    }

    (indices, values)
}

/// Build a table of Ballot numbers B_pq from p=q=0 to p=q=b.
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

/// Compute the block types for all blocks in x and compute the tables for the
/// blocks we observe.
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

pub struct Optimal<'a> {
    x: &'a [usize],
    block_size: BlockSize,
    sparse: Sparse,
    block_types: Vec<usize>,
    block_tables: Vec<Option<TabulatedQuery>>,
}

impl<'a> Optimal<'a> {
    pub fn new(x: &'a [usize]) -> Self {
        let n = x.len();
        let BlockSize(b) = block_size(n);

        // adjust block size; log(n) is too much so change it to a quarter.
        let b = std::cmp::max(4, b / 4); // I don't want too small blocks, so minimum is 4
        let block_size = BlockSize(b);

        let (_, reduced_vals) = reduce_array(x, block_size);
        let (block_types, block_tables) = tabulate_blocks(x, b);

        Optimal {
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

        // Get the table for this block by looking up the block type and then the
        // table from the block type.
        let btype = self.block_types[block_index];
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

impl<'a> RMQ for Optimal<'a> {
    fn rmq(&self, i: usize, j: usize) -> Option<usize> {
        let BlockSize(bs) = self.block_size;
        // The block indices are not the same for the small tables and the
        // sparse table. For the sparse table we have to round up for i, but
        // to get the block i is in, we need to round down.
        let bi = BlockIdx(i / bs);
        let (sparse_bi, ii) = round_up(i, BlockSize(bs));
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
