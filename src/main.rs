use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

#[inline(always)]
fn flog2(v: usize) -> usize {
    v.ilog2() as usize
}

trait Rmq {
    fn rmq(&self, i: usize, j: usize) -> usize;
}

struct Sparse {
    lgn: usize,
    table: Vec<usize>, // convert this to a slice
}

impl Sparse {
    fn new(t: &[usize]) -> Self {
        let n = t.len();
        let lgn = flog2(n);

        let mut table: Vec<usize> = vec![0; (lgn + 1) * n];

        for j in 0..n {
            table[j * (lgn + 1)] = t[j];
        }

        // Build the sparse table
        for i in 1..=lgn {
            for j in 0..=(n - (1 << i)) {
                table[j * (lgn + 1) + i] = std::cmp::min(
                    table[j * (lgn + 1) + i - 1],
                    table[(j + (1 << (i - 1))) * (lgn + 1) + i - 1],
                );
            }
        }

        Sparse { lgn, table }
    }
}

impl Rmq for Sparse {
    fn rmq(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < j);

        let k = flog2(j - i + 1);
        std::cmp::min(
            self.table[i * (self.lgn + 1) + k],
            self.table[(j + 1 - (1 << k)) * (self.lgn + 1) + k],
        )
    }
}

fn create_random_ary(ary_size: usize, max_ary: usize) -> Vec<usize> {
    let mut rng = StdRng::from_entropy();
    (0..ary_size).map(|_| rng.gen_range(0..=max_ary)).collect()
}

fn main() {
    let ary_size = 10;
    let max_ary = 100;
    let ary = create_random_ary(ary_size, max_ary);
    let sparse = Sparse::new(&ary);

    for i in 0..ary_size {
        for j in i + 1..ary_size {
            let res = sparse.rmq(i, j);
            let min_naive = (i..=j).map(|k| ary[k]).min().unwrap();
            assert_eq!(res, min_naive)
        }
    }

    println!("ALL OK!");
}

// mod tests {
//     use crate::*;

//     #[test]
//     fn test_correct() {
//         let ary_size = 100;
//         let max_ary = 1000;
//         let ary = create_random_ary(ary_size, max_ary);
//         let sparse = Sparse::new(&ary);

//         for i in 0..ary_size {
//             for j in i+1..ary_size {
//                 let res = sparse.rmq(i, j);
//                 let min_naive = (i+1..=j).map(|k| ary[k]).min().unwrap();
//                 assert_eq!(ary[res], min_naive)
//             }
//         }
//     }
// }
