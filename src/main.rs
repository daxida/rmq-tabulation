mod rmq;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rmq::{Rmq, Sparse};

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
