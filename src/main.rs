use lib::{Cartesian, Sparse, Rmq};

fn compare_methods(ary: &[usize]) {
    let ary_size = ary.len();
    let sparse = Sparse::new(ary);
    let opt = Cartesian::new(ary);

    for i in 0..ary_size {
        for j in i + 1..ary_size {
            let res = sparse.rmq(i, j).unwrap();
            let min_naive = (i..j).map(|k| ary[k]).min().unwrap();
            let res_opt = opt.rmq(i, j).unwrap_or(0);
            assert_eq!(res, min_naive);
            // assert_eq!(res_sparse_opt, min_naive);
            assert_eq!(res_opt, min_naive);
            // println!("{} {} {} {}", i, j, res, res_opt);
        }
    }
}

fn main() {
    let ary: Vec<usize> = vec![
        251, // 0  <- min idx 0
        339, // 1
        458, // 2
        261, // 3  BLOCK
        389, // 4
        375, // 5  <- min idx 5
        908, // 6
        888, // 7  BLOCK
        697, // 8
        676, // 9
        169, // 10 <- min idx 10
        607, // 11 BLOCK
        912, // 12
        160, // 13
        281, // 14
        44,  // 15 BLOCK <- min idx 15
    ];

    compare_methods(&ary);

    println!("ALL OK!");
}

#[cfg(test)]
mod tests {
    use super::*;

    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

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
}
