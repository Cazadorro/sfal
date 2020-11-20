/// Functions that approximate the erfinv(x), the inverse "Error Function".

use super::consts;
use std::f64;

///https://en.wikipedia.org/wiki/Error_function#Inverse_functions
pub fn maclaurin_series<const MAX_N: usize>(x:f64) -> f64{
    let mut ck_array = [0.0; MAX_N];
    ck_array[0] = 1.0;
    let mut product = 1.0;
    for k in 1..MAX_N{
        let mut ck = 0.0;
        for m in 0..=(k-1){
            ck += (ck_array[m as usize]*ck_array[(k - 1 - m) as usize])/((m+1)*(2*m + 1)) as f64;
        }
        product *= (f64::consts::PI / 4.0);
        ck_array[k] = ck * product;
    }
    let mut sum = x;
    let mut x_prod = x;
    for i in 1..MAX_N{
        x_prod *= x * x;
        sum += ck_array[i] * x_prod;
    }
    (f64::consts::PI.sqrt()/2.0) * sum
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn sergei_pade(x:f64)-> f64{
    let x_negative = x < 0.0;
    let x = if x_negative { -x } else { x };

    let a = (8.0 * (f64::consts::PI - 3.0))/ ( 3.0*f64::consts::PI*(3.0 - f64::consts::PI));
    let ln1_x2 = f64::ln(1.0 - x*x);
    let inner_root = (2.0/(f64::consts::PI*a)) + (ln1_x2/2.0);
    let result = inner_root * inner_root;
    let result = result - (ln1_x2/a);
    let result = result.sqrt();
    let result = result - inner_root;
    let result = x.signum() * result.sqrt();

    if x_negative {
        -result
    } else {
        result
    }
}

///http://www.mimirgames.com/articles/programming/approximations-of-the-inverse-error-function/#
pub fn newton_refinement(x : f64, refinement_count : u64, erf_fn : impl Fn(f64) -> f64) -> f64{

    let mut r = 0.0;
    let a = [0.886226899, -1.645349621, 0.914624893, -0.140543331];
    let b = [1.0, -2.118377725, 1.442710462, -0.329097515, 0.012229801];
    let c = [-1.970840454, -1.62490649, 3.429567803, 1.641345311];
    let d = [1.0, 3.543889200, 1.637067800];

    let z = x.signum() * x;

    if z <= 0.7{
        let x2 = z * z;
        r = z * (((a[3] * x2 + a[2]) * x2 + a[1]) * x2 + a[0]);
        r /= (((b[4] * x2 + b[3]) * x2 + b[2]) * x2 + b[1])* x2 + b[0];
    } else {
        let y = f64::sqrt( -f64::ln((1.0 - z)/2.0));
        r = (((c[3] * y + c[2]) * y + c[1]) * y + c[0]);
        r /= ((d[2] * y + d[1]) * y + d[0]);
    }

    r = r * x.signum();
    let z = z * x.signum();

    for i in  0..=refinement_count{
        r -= (erf_fn(r) - z)/(2.0/f64::consts::PI.sqrt() *f64::exp(-r * r));
    }
    r
}


///http://www.mimirgames.com/articles/programming/approximations-of-the-inverse-error-function/#
///Numerical Recipes Third Edition p265
pub fn halley_refinement(x : f64, refinement_count : u64, erfc_fn : impl Fn(f64) -> f64) -> f64{
    let mut pp= 0.0;
    let mut t = 0.0;
    let mut r = 0.0;
    let mut err = 0.0;

    if x < 1.0{
        pp = x;
    } else {
        pp = 2.0 - x;
    }
    t = f64::sqrt(-2.0 * f64::ln(pp/2.0));
    r = -f64::consts::FRAC_1_SQRT_2 * ((2.30753 + t * 0.27061)/(1.0 + t * (0.99229 + t * 0.04481)) - t);

    for i in  0..=refinement_count {
        err = erfc_fn(r) - pp;
        r += err/(f64::consts::FRAC_2_SQRT_PI * f64::exp(-r * r) - r * err);
    }
    if x > 1.0{
        r = -r;
    }
    r
}


///from "Approximating the erfinv function" by Mike Giles
pub fn mike_giles(x : f64) -> f64{
    let mut p = 0.0;
    let mut w = -f64::ln((1.0 - x)*(1.0 + x));
    if w < 0.5{
        w = w - 2.5;
        p = 2.81022636e-08;
        p = 3.43273939e-07 + p * w;
        p = -3.5233877e-06 + p * w;
        p = -4.39150654e-06 + p * w;
        p = 0.00021858087 + p * w;
        p = -0.00125372503 + p * w;
        p = -0.00417768164 + p * w;
        p = 0.246640727 + p * w;
        p = 1.50140941 + p * w;
    } else {
        w = w.sqrt() - 3.0;
        p = -0.000200214257;
        p = 0.000100950558 + p * w;
        p = 0.00134934322 + p * w;
        p = -0.00367342844 + p * w;
        p = 0.00573950773 + p * w;
        p = -0.0076224613 + p * w;
        p = 0.00943887047 + p * w;
        p = 1.00167406 + p * w;
        p = 2.83297682 + p * w;
    }
    return p * x;
}