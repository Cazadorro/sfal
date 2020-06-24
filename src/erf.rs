/// Functions that approximate the erf(x), the "Error Function".

use super::consts;
use std::f64;

///https://en.wikipedia.org/wiki/Error_function
fn convergent_series(x:f64, max_n:u64) -> f64{
    let mut total_sum = 0.0;
    let mut n_coefficient = x;
    for n in 1..max_n{

        n_coefficient *= -(x*x);
        n_coefficient /= n as f64;

        let prev_sum = total_sum;
        total_sum += n_coefficient * (1.0/(2*n + 1) as f64);
        //comparison against itself basically, shouldn't run into pitfalls with float equality.
        if total_sum == prev_sum{
            break;
        }
    }
    2.0/(std::f64::consts::PI.sqrt()) * total_sum
}

///https://en.wikipedia.org/wiki/Error_function#B%C3%BCrmann_series
fn burmann_short(x:f64, max_n:u64) -> f64{
3.0
}

