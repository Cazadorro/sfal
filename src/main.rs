
#![feature(const_generics)]
#![feature(test)]
// #![feature(fn_traits)]
// #![feature(unboxed_closures)]
/// calculates erf using series expansion shown here:
/// http://people.math.sfu.ca/~cbm/aands/page_297.htm
// fn erf_series_0()
use std::f64;
use std::ops;
use num::complex::{Complex32,Complex64};
mod consts;
mod cheb;
mod ei;
mod e1;

//
// extern crate test;

pub trait MiscMath<ReturnType = f64, IntegralType=u64> {
    fn powi_sign(self) -> ReturnType;
    fn factorial(self) -> Self;
    fn offset_factorial(self, offset:Self) -> Self;
}

//https://en.wikipedia.org/wiki/Exponential_integral#Approximations
fn E1_Swammee_Ohija(x:f64) -> f64{
    let A = f64::ln((0.56146/x + 0.65) * (1.0 + x));
    let B = x.powi(4)* f64::exp(7.7*x) * (2.0+x).powf(3.7);
    (A.powf(-7.7) + B).powf(-0.13)
}

fn E1_Allen_Hastings(x:f64) -> f64{
    let a = [-0.57722, 0.99999, -0.24991, 0.05519, -0.00976, 0.00108];
    let b = [0.26777, 8.63476, 18.05902, 8.57333];
    let c = [3.95850, 21.09965, 25.63296, 9.57332];
    let x3 = [1.0, x, x.powi(2), x.powi(3)];
    let x5 = [1.0, x, x.powi(2), x.powi(3), x.powi(4), x.powi(5)];
    if x <= 1.0{
        let aTx5 = a.iter().zip(&x5).map(|(ai, xi)| (*ai)*(*xi)).sum::<f64>();
        -x.ln() + aTx5
    }else{
        let bTx3 = (b.iter().zip(&x3).map(|(bi, xi)|  (*bi)*(*xi)).sum::<f64>());
        let cTx3 = (c.iter().zip(&x3).map(|(ci, xi)| (*ci)*(*xi)).sum::<f64>());
        (f64::exp(-x)/x) * (bTx3/cTx3)
    }
}

fn E1_continued_fraction(x:f64, iterations:u64) -> f64{
    let mut fraction = 1.0;
    for i in (1..=iterations).rev(){
       fraction = x + (i as f64/(1.0+(i as f64/fraction)));
    }
    f64::exp(-x) / fraction
}

fn E1_Barry_et_al(x:f64) -> f64{
    let q = (20.0/47.0) * x.powf((31.0_f64/26.0).sqrt());
    let G = (-consts::EULER_MASCHERONI).exp();
    let b = ((2.0*(1.0 - G))/(G*(2.0 - G))).sqrt();
    let h_inf = ((1.0 - G) * (G*G - 6.0*G + 12.0))/(3.0*G * (2.0 - G).powi(2)*b);
    let h = (1.0/(1.0 + x * x.sqrt())) + (h_inf*q)/(1.0 +q);
    ((-x).exp()/(G + (1.0 - G) * f64::exp(-x/(1.0-G)))) * f64::ln((1.0 + G/x - (1.0 - G)/(h + b*x).powi(2)).abs())
}

fn E1_Barry_et_al2(x:Complex64) -> Complex64{
    let q = (20.0/47.0) * x.powf((31.0_f64/26.0).sqrt());
    let G = (-consts::EULER_MASCHERONI).exp();
    let b = ((2.0*(1.0 - G))/(G*(2.0 - G))).sqrt();
    let h_inf = ((1.0 - G) * (G*G - 6.0*G + 12.0))/(3.0*G * (2.0 - G).powi(2)*b);
    let h = (1.0/(1.0 + x * x.sqrt())) + (h_inf*q)/(1.0 +q);
    ((-x).exp()/(G + (1.0 - G) * Complex64::exp(-x/(1.0-G)))) * Complex64::ln(1.0 + G/x - (1.0 - G)/(h + b*x).powi(2))
}

fn calc_u(x: f64) -> f64 {
    return (x - 1.0) / 2.0;
}

fn calc_x(u: f64) -> f64 {
    return 2.0 * u + 1.0;
}

fn chebyshev_node(i: u64, n: u64) -> f64 {
    return (((2.0 * i as f64 - 1.0) * f64::consts::PI) / (2.0 * n as f64)).cos();
}

fn fx(x: f64) -> f64 {
    return (1.0 / 3.0) * x.powi(3) + 2.0 * x.powi(2) + x - 10.0;
}

fn erf_0(x: f32) -> f32 {
    let x = x.abs();
    let a1 = 0.278393_f32;
    let a2 = 0.230389_f32;
    let a3 = 0.000972_f32;
    let a4 = 0.078108_f32;
    let x2 = x * x;
    let x3 = x2 * x;
    let x4 = x3 * x;

    let denominator = 1_f32 + a1 * x + a2 * x2 + a3 * x3 + a4 * x4;
    let denominator = denominator.powi(4);

    return 1_f32 - (1_f32 / denominator);
}

#[cfg(test)]
mod tests{
    use super::*;
    use test::Bencher;

    #[bench]
    fn bench_stuff(b: &mut Bencher){
        b.iter(|| Ei_series(24.345, 28));
    }
}
fn main() {
    let mut nodes = [0.0; 5];
    (1..=5_u64).for_each(|i| nodes[(i - 1) as usize] = chebyshev_node(i, 5));
    println!("nodes");
    let nodes = nodes;
    (0..5_usize).for_each(|i| println!("{}", nodes[i]));

    println!("converted");
    let mut xs = [0.0; 5];
    (0..5_usize).for_each(|i| xs[i] = calc_x(nodes[i]));
    (0..5_usize).for_each(|i| println!("{}", xs[i]));

    println!("pure");
    let mut ys0 = [0.0; 5];
    (0..5_usize).for_each(|i| ys0[i] = fx(nodes[i]));
    (0..5_usize).for_each(|i| println!("{}", ys0[i]));

    println!("converted");
    let mut ys1 = [0.0; 5];
    (0..5_usize).for_each(|i| ys1[i] = fx(xs[i]));
    (0..5_usize).for_each(|i| println!("{}", ys1[i]));

    println!("{}", erf_0(3.0));
    println!("{}", erf_0(3.0));
    println!("Hello, world!");

    let test = cheb::gen_chebyshev_approximator(-1.0,3.0,4_u64,(|x| (1.0 / 3.0) * x.powi(3) + 2.0 * x.powi(2) + x - 10.0));
    // println!("{:#?}", test);
    println!("{}", test(0.0));

    let ei_cheb = cheb::gen_chebyshev_approximator(0.2,2.0,32_u64,(|x| ei::convergent_ramanujan(x, 28)));

    // E1_Swammee_Ohija
    // E1_Allen_Hastings
    // E1_continued_fraction
    // E1_Barry_et_al


    const APPROX_LEN :usize = 20;
    let mut eix_array = [0.0; APPROX_LEN];
    let mut ei_array = [0.0; APPROX_LEN];
    let mut ei_cheb_array = [0.0; APPROX_LEN];
    let mut ei_swamee_array = [0.0; APPROX_LEN];
    let mut ei_allen_array = [0.0; APPROX_LEN];
    let mut ei_continued_4 = [0.0; APPROX_LEN];
    let mut ei_continued_16 = [0.0; APPROX_LEN];
    let mut ei_barry_array = [0.0; APPROX_LEN];
    let mut ei_barry_array2 = [0.0; APPROX_LEN];
    let mut e1_array = [0.0; APPROX_LEN];
    let mut e1_continued_4 = [0.0; APPROX_LEN];
    let mut e1_continued_16 = [0.0; APPROX_LEN];
    let mut ei_large_small = [0.0; APPROX_LEN];
    let mut ei_large = [0.0; APPROX_LEN];
    let mut ei_small = [0.0; APPROX_LEN];
    let mut ei_array_2 = [0.0; APPROX_LEN];
    let mut ei_array_opt = [0.0; APPROX_LEN];
    let mut e1_array_opt = [0.0; APPROX_LEN];
    let mut nei_array_opt = [0.0; APPROX_LEN];
    let mut ei_divergent_test = [0.0; APPROX_LEN];
    let mut ei_divergent_test2 = [0.0; APPROX_LEN];


    for i in 0..APPROX_LEN{

        eix_array[i] = i as f64 * 2.0/(ei_array.len() as f64);
        let x = eix_array[i];
        ei_array[i] = ei::convergent_ramanujan(x, 28);
        ei_cheb_array[i] = ei_cheb(x);
        // ei_swamee_array[i] = -E1_Swammee_Ohija(-x);
        // ei_allen_array[i] = -E1_Allen_Hastings(-x);
        // ei_continued_4[i] = -E1_continued_fraction(-x, 4);
        // ei_continued_16[i] = -E1_continued_fraction(-x, 16);
        // ei_barry_array[i] = -E1_Barry_et_al(-x);
        // ei_barry_array2[i] = -E1_Barry_et_al2( Complex64::from(-x)).re;
        // e1_array[i] = E1_Abramowitz_Stegun_series(x, 28);
        // e1_continued_4[i] = E1_continued_fraction(x, 4);
        // e1_continued_16[i] = E1_continued_fraction(x, 16);
        let x2 = x*10.0;
        // ei_large_small[i] = Ei_large_small(x2, 30);
        // ei_large[i] = Ei_large(x2, 100);
        // ei_small[i] = Ei_small(x2, 40);
        ei_array_2[i] = ei::convergent_ramanujan(x2, 64);
        ei_array_opt[i] = ei::convergent_ramanujan(x2, 64);
        e1_array_opt[i] = ei::convergent_ramanujan(x2, 64);
        nei_array_opt[i] =  ei::convergent_ramanujan(-x, 64);
        // ei_divergent_test[i] = ei::divergent_series(-30.0, (((i+1))+25) as u64);
        ei_divergent_test[i] = -ei::continued_fraction(-x2, 1023 as u64);
        ei_divergent_test2[i] = e1::continued_fraction(19.0, ((i+1)+40) as u64);
    }

    println!("{:#?}", eix_array);
    println!("{:#?}", ei_array);
    println!("{:#?}", ei_cheb_array);
    // println!("{:#?}", ei_swamee_array);
    // println!("{:#?}", ei_allen_array);
    // println!("{:#?}", ei_continued_4);
    // println!("{:#?}", ei_continued_16);
    println!("{:#?}", ei_barry_array);
    println!("{:#?}", ei_barry_array2);
    println!("{:#?}", e1_array);
    println!("{:#?}", e1_continued_4);
    println!("{:#?}", e1_continued_16);
    println!("{:#?}", ei::convergent_ramanujan(2.0, 28));
    println!("{:#?}", ei::convergent_ramanujan(5.0, 28));
    println!("{:#?}", -E1_Barry_et_al2( Complex64::from(-0.1)));
    println!("{:#?}", ei_large_small);
    println!("{:#?}", ei_small);
    println!("{:#?}", ei_large);
    println!("{:#?}", ei_array_2);
    println!("{:#?}", ei_array_opt);
    println!("{:#?}", e1_array_opt);
    println!("{:#?}", e1_array);
    println!("{:#?}", nei_array_opt);
    println!("{:#?}", ei_divergent_test);
    println!("{:#?}", ei_divergent_test2);
}


