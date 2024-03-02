mod range;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use range::{Rmq, SparseT, RMQ, Optimal};

fn create_random_ary(ary_size: usize, max_ary: usize) -> Vec<usize> {
    let mut rng = StdRng::from_entropy();
    (0..ary_size).map(|_| rng.gen_range(0..=max_ary)).collect()
}

fn main() {
    let ary_size = 50;
    let max_ary = 1000;
    let ary = create_random_ary(ary_size, max_ary);
    let sparse = SparseT::new(&ary);
    let opt = Optimal::new(&ary);

    for i in 0..ary_size {
        for j in i + 1..ary_size {
            let res = sparse.rmq(i, j);
            let min_naive = (i..=j).map(|k| ary[k]).min().unwrap();
            let res_opt = ary[opt.rmq(i, j+1).unwrap_or(0)];
            assert_eq!(res, min_naive);
            assert_eq!(res_opt, min_naive);
        }
    }

    println!("ALL OK!");
}