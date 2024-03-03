use lib::{Optimal, Rmq, SparseT, RMQ};

fn compare_methods(ary: &Vec<usize>) {
    let ary_size = ary.len();
    let sparse = SparseT::new(&ary);
    let opt = Optimal::new(&ary);

    for i in 0..ary_size {
        for j in i + 1..ary_size {
            let res = sparse.rmq(i, j);
            let min_naive = (i..=j).map(|k| ary[k]).min().unwrap();
            let res_opt = opt.rmq(i, j + 1).unwrap_or(0);
            assert_eq!(res, min_naive);
            assert_eq!(res_opt, min_naive);
        }
    }
}

fn main() {
    let ary: Vec<usize> = vec![
        339,
        458,
        261,
        251, // <- min idx 3
        389, //
        908,
        375, // <- min idx 6
        888,
        697,
        676, //
        169, // <- min idx 10
        607,
        912,
        160,
        281, //
        188,
        890,
        681,
        51,  // <- min idx 18
        947, //
        162, // <- min idx 20
        891,
        705,
        755,
        877,
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
        let ary_size = 500;
        let max_ary = 10;
        let ary = create_random_ary(ary_size, max_ary);
        compare_methods(&ary)
    }
}
