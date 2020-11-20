/// Functions that approximate the erf(x), the "Error Function".

use super::erfc;
use super::consts;
use std::f64;

///https://en.wikipedia.org/wiki/Error_function
pub fn convergent_series(x:f64, max_n:u64) -> f64{
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
pub fn burmann_short(x:f64) -> f64{
    let x2 = x*x;
    let sum = (f64::consts::PI.sqrt()/2.0) + (31.0/200.0)*f64::exp(-x2) - (341.0/8000.0)*f64::exp(-2.0*x2);
    (2.0/f64::consts::PI.sqrt()) * x.signum() * f64::sqrt(1.0 - f64::exp(-(x2))) * sum
}

///https://en.wikipedia.org/wiki/Error_function#Asymptotic_expansion
pub fn divergent_series(x:f64, max_n:u64) -> f64{
    1.0 - erfc::divergent_series(x, max_n)
}

///https://en.wikipedia.org/wiki/Error_function#Continued_fraction_expansion
pub fn continued_fraction(x:f64, max_n:u64) -> f64{
   1.0 - erfc::continued_fraction(x, max_n)
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn abramowitz_0(x: f64) -> f64 {
    1.0 - erfc::abramowitz_0(x)
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn abramowitz_1(x: f64) -> f64 {
    1.0 - erfc::abramowitz_1(x)
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn abramowitz_2(x: f64) -> f64 {
    1.0 - erfc::abramowitz_2(x)
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn abramowitz_3(x: f64) -> f64 {
    1.0 - erfc::abramowitz_3(x)
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn karagiannidis(x: f64) -> f64{
   1.0 - erfc::karagiannidis(x)
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn sergei_pade(x:f64)-> f64{
    let x_negative = x < 0.0;
    let x = if x_negative { -x } else { x };

    let a = (8.0 * (f64::consts::PI - 3.0))/ ( 3.0*f64::consts::PI*(3.0 - f64::consts::PI));
    let ax2 = a *(x*x);
    let result = -(x*x) * ((4.0/f64::consts::PI) + ax2)/(1.0 + ax2);
    let result = 1.0 - result.exp();
    let result = x.signum()*result.sqrt();
    if x_negative {
        -result
    } else {
        result
    }
}


///https://en.wikipedia.org/wiki/Error_function#Polynomial
pub fn polynomial_9th(x:f64) -> f64{
    let ai_array = [-1.26551223, 1.00002368, 0.37409196, 0.09678418, -0.18628806, 0.27886807, -1.13520398, 1.48851587, -0.82215223, 0.17087277];
    let t = 1.0 / (1.0 + (x.abs()/2.0));
    let mut sum = x*x + ai_array[0];
    let mut t_prod = 1.0;
    for i in 1..ai_array.len(){
        t_prod *= t;
        sum += ai_array[i]*t_prod;
    }
    let tau = t*f64::exp(sum);
    if x >= 0.0{
        1.0 - tau
    }else{
        tau - 1.0
    }
}

///Numerical Recipes third edition page 265.
pub fn recipes_cheb(x : f64) -> f64{
    1.0 - erfc::recipes_cheb(x)
}