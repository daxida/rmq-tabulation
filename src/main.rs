use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;


#[inline(always)]
fn flog2(v: usize) -> usize {
    v.ilog2() as usize
}

trait Rmq {
    fn rmq(&self, i: usize, j: usize) -> usize;
}

struct Sparse<'a> {
    t: &'a [usize],
    n: usize,
    table: Vec<usize>, // convert this to a slice
}

#[inline(always)]
fn min_by_t(t: &[usize], l: usize, r: usize) -> usize {
    if t[l] < t[r] {
        l
    } else {
        r
    }
}

impl<'a> Sparse<'a> {
    fn new(t: &'a [usize]) -> Self {
        let n = t.len();
        let lgn = flog2(n);

        let mut table = vec![0; n * lgn];
        for i in 0..n {
            table[i * lgn] = i;
        }

        let mut j = 1;
        while (1 << j) <= n {
            for i in 0..=n - (1 << j) {
                let idx_l = i * lgn + j;
                let idx_r = (i + (1 << (j - 1))) * lgn + j - 1;
                let l = table[idx_l - 1];
                let r = table[idx_r];
                table[idx_l] = min_by_t(t, l, r);
            }
            j += 1;
        }

        Sparse { t, n, table }
    }
}

impl <'a> Rmq for Sparse<'a> {
    fn rmq(&self, mut i: usize, j: usize) -> usize {
        debug_assert!(i < j);
        i += 1;
    
        if i < j {
            let lgn = flog2(self.n);
            let k = flog2(j - i + 1);
    
            // Calculate indices a and b for the two halves of the range.
            let idx_l = i * lgn + k;
            let idx_r = (j - (1 << k) + 1) * lgn + k;
            let l = self.table[idx_l];
            let r = self.table[idx_r];
    
            min_by_t(self.t, l, r)
        } else {
            // i == j since i <= j
            i
        }  
    }
}

fn create_random_ary(ary_size: usize, max_ary: usize) -> Vec<usize> {
    let mut rng = StdRng::from_entropy();
    (0..ary_size).map(|_| rng.gen_range(0..=max_ary)).collect()
}

fn main() {
    let ary_size = 100;
    let max_ary = 1000;
    let ary = create_random_ary(ary_size, max_ary);
    let sparse = Sparse::new(&ary);

    for i in 0..ary_size {
        for j in i+1..ary_size {
            let res = sparse.rmq(i, j);
            let min_naive = (i+1..=j).map(|k| ary[k]).min().unwrap();
            assert_eq!(ary[res], min_naive)
        }
    }
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