/// Functions that approximate the erfc(x), the "Complementary Error Function".

use super::erf;
use super::consts;
use std::f64;

///https://en.wikipedia.org/wiki/Error_function#Asymptotic_expansion
pub fn divergent_series(x: f64, max_n: u64) -> f64 {
    let mut prod = 1.0;
    let mut sum = 1.0;
    for n in 1..max_n {
        if n & 1 == 1 {
            prod *= n as f64;
        }
        prod /= (2.0 * x * x);
    }
    (f64::exp(-(x * x)) / (x * f64::consts::PI.sqrt())) * sum
}

///https://en.wikipedia.org/wiki/Error_function#Continued_fraction_expansion
pub fn continued_fraction(x: f64, max_n: u64) -> f64 {
    let mut fraction = 1.0;
    for i in (1..=max_n).rev() {
        let ai = i as f64 / 2.0;
        if i & 1 == 1 {
            fraction = (x * x) + (ai as f64 / fraction);
        } else {
            fraction = 1.0 + (ai as f64 / fraction);
        }
    }
    (x / f64::consts::PI.sqrt()) * f64::exp(-(x * x)) / fraction
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn abramowitz_0(x: f64) -> f64 {
    let x_negative = x < 0.0;
    let x = if x_negative { -x } else { x };
    let ai_array = [0.278393, 0.230389, 0.000972, 0.078108];
    let mut denominator_sum = 1.0_f64;
    let mut x_prod = 1.0;
    for i in 0..ai_array.len() {
        x_prod *= x;
        denominator_sum += ai_array[i] * x_prod;
    }
    let result = (1.0 / denominator_sum.powi(4));
    if x_negative {
        -result
    } else {
        result
    }
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn abramowitz_1(x: f64) -> f64 {
    let x_negative = x < 0.0;
    let x = if x_negative { -x } else { x };
    let p = 0.47047;
    let t = 1.0/(1.0 + p*x);
    let ai_array = [0.3480242,-0.0958798,0.7478556];
    let mut sum = 1.0;
    let mut t_prod = 1.0;
    for i in 0..ai_array.len() {
        t_prod *= t;
        sum += ai_array[i] * t_prod;
    }
    let result = (sum*f64::exp(-(x*x)));
    if x_negative {
        -result
    } else {
        result
    }
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn abramowitz_2(x: f64) -> f64 {
    let x_negative = x < 0.0;
    let x = if x_negative { -x } else { x };
    let ai_array = [0.0705230784, 0.0422820123,  0.0092705272, 0.0001520143, 0.0002765672, 0.0000430638];
    let mut denominator_sum = 1.0_f64;
    let mut x_prod = 1.0;
    for i in 0..ai_array.len() {
        x_prod *= x;
        denominator_sum += ai_array[i] * x_prod;
    }

    let result = (1.0 / denominator_sum.powi(16));
    if x_negative {
        -result
    } else {
        result
    }
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn abramowitz_3(x: f64) -> f64 {
    let x_negative = x < 0.0;
    let x = if x_negative { -x } else { x };
    let p = 0.3275911;
    let t = 1.0/(1.0 + p*x);
    let ai_array = [0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429];
    let mut sum = 1.0;
    let mut t_prod = 1.0;
    for i in 0..ai_array.len() {
        t_prod *= t;
        sum += ai_array[i] * t_prod;
    }
    let result = (sum*f64::exp(-(x*x)));
    if x_negative {
        -result
    } else {
        result
    }
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn karagiannidis(x: f64) -> f64{
    let x_negative = x < 0.0;
    let x = if x_negative { -x } else { x };


    let a = 1.98;
    let b = 1.135;

    let result = ((1.0 - f64::exp(-a*x))*f64::exp(-(x*x)))/(b * f64::consts::PI.sqrt() * x);
    if x_negative {
        -result
    } else {
        result
    }
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn sergei_pade(x:f64)-> f64{
    1.0 - erf::sergei_pade(x)
}

pub fn polynomial_9th(x:f64) -> f64{
   1.0 - erf::polynomial_9th(x)
}

///Numerical Recipes third edition page 265.
pub fn recipes_cheb(x : f64) -> f64{
    let cheb_coef_array = [-1.3026537197817094, 6.4196979235649026e-1,
        1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
        3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
        -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
        6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
        9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
        -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17];

    let mut d=0.0;
    let mut dd=0.0;
    let x_negative = x < 0.0;
    let x = if x_negative{ -x }else{ x };

    let t = 2.0/(2.0+x);
    let ty = 4.0*t - 2.0;
    for j in (0..cheb_coef_array.len()).rev() {
        let tmp = d;
        d = ty*d - dd + cheb_coef_array[j as usize];
        dd = tmp;
    }
    let result = t*f64::exp(-x*x + 0.5*(cheb_coef_array[0] + ty*d) - dd);
    if x_negative{
        -result
    }else{
        result
    }
}
