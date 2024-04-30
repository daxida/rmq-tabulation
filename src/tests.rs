use super::{Rmq, Sparse, Tabulation};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn compare_methods(ary: &[usize]) {
    let ary_size = ary.len();
    let sparse = Sparse::new(ary);
    let opt = Tabulation::new(ary);

    for i in 0..ary_size {
        for j in i + 1..ary_size {
            let res = sparse.rmq(i, j).unwrap();
            let min_naive = (i..j).map(|k| ary[k]).min().unwrap();
            let res_opt = opt.rmq(i, j).unwrap_or(0);
            assert_eq!(res, min_naive);
            assert_eq!(res_opt, min_naive);
        }
    }
}

fn create_random_ary(ary_size: usize, max_ary: usize) -> Vec<usize> {
    let mut rng = StdRng::from_entropy();
    (0..ary_size).map(|_| rng.gen_range(0..=max_ary)).collect()
}

#[test]
fn test_random_one() {
    let ary_size = 20;
    let max_ary = 1000;
    let ary = create_random_ary(ary_size, max_ary);
    compare_methods(&ary)
}

#[test]
fn test_random_two() {
    let ary_size = 200;
    let max_ary = 1000000;
    let ary = create_random_ary(ary_size, max_ary);
    compare_methods(&ary)
}

#[test]
fn test_random_three() {
    let ary_size = 250;
    let max_ary = 10;
    let ary = create_random_ary(ary_size, max_ary);
    compare_methods(&ary)
}

#[test]
fn test_random_four() {
    let ary_size = 1 << 4;
    let max_ary = 10;
    let ary = create_random_ary(ary_size, max_ary);
    compare_methods(&ary)
}

#[test]
fn test_random_five() {
    let ary_size = 1 << 5;
    let max_ary = 5;
    let ary = create_random_ary(ary_size, max_ary);
    compare_methods(&ary)
}
