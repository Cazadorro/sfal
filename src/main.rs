

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
mod errors;
mod cheb;
mod burmann;
mod ei;
mod e1;
mod erf;
mod erfc;
mod erfcinv;
mod erfinv;

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

// #[cfg(test)]
// mod tests{
//     use super::*;
//     use test::Bencher;
//
//     #[bench]
//     fn bench_stuff(b: &mut Bencher){
//         b.iter(|| Ei_series(24.345, 28));
//     }
// }

fn create_cos(node_count:u64, divider:f64) -> impl Fn(f64) -> f64{
    let eighth_circle = f64::consts::PI/4.0;
    let quarter_circle = f64::consts::PI/2.0;
    let small_angle = 0.0000;
    let cos_gen = cheb::gen_chebyshev_approximator(small_angle, 1.0/divider, node_count, (|x| f64::cos(x*f64::consts::PI)));
    let approximator = move |x: f64| -> f64{
        let x = if x < 0.0 { -x }else{ x };
        if x < small_angle{
            1.0 - (x*x)/2.0
        }else{
            cos_gen(x/f64::consts::PI)
        }
        // let x = if x < 0.0 { -x }else{ x };
        // let x = x % 2.0*f64::consts::PI;
        // let x = if x > f64::consts::PI { 2.0*f64::consts::PI  - x } else { x };
        // let invert_output = x > quarter_circle;
        // if x > eighth_circle {
        //     let x = 2.0* eighth_circle - x;
        //     let y = 1.0 - cos_gen(x);
        //     if invert_output{
        //         -y
        //     }else{
        //         y
        //     }
        //
        // }else {
        //     let y = cos_gen(x);
        //     if invert_output{
        //         -y
        //     }else{
        //         y
        //     }
        // }
    };
    approximator
}

fn bhaskara(x:f64)->f64{
    ((f64::consts::PI.powi(2) - 4.0 * x.powi(2))/(f64::consts::PI.powi(2) + x.powi(2)))
}
fn test(x:bool, f: fn(bool) -> bool) -> bool{
    f(x)
}

fn create_cos2(node_count:u64, divider:f64) -> impl Fn(f64) -> f64{
    let eighth_circle = f64::consts::PI/4.0;
    let quarter_circle = f64::consts::PI/2.0;
    let small_angle = 0.1;
    let cos_gen = cheb::gen_chebyshev_approximator(small_angle, 1.0/divider, node_count, (|x|  (f64::cos(x) - bhaskara(x))));
    let approximator = move |x: f64| -> f64{
        let x = if x < 0.0 { -x }else{ x };
        if x < small_angle{
            1.0 - (x*x)/2.0
        }else{
            bhaskara(x) + cos_gen(x)
        }

        // cos_gen(x/f64::consts::PI)
        // let x = if x < 0.0 { -x }else{ x };
        // let x = x % 2.0*f64::consts::PI;
        // let x = if x > f64::consts::PI { 2.0*f64::consts::PI  - x } else { x };
        // let invert_output = x > quarter_circle;
        // if x > eighth_circle {
        //     let x = 2.0* eighth_circle - x;
        //     let y = 1.0 - cos_gen(x);
        //     if invert_output{
        //         -y
        //     }else{
        //         y
        //     }
        //
        // }else {
        //     let y = cos_gen(x);
        //     if invert_output{
        //         -y
        //     }else{
        //         y
        //     }
        // }
    };
    approximator
}
fn main() {
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
    let mut e1_swamee = [0.0; APPROX_LEN];
    let mut ei_swamee = [0.0; APPROX_LEN];
    let mut ei_swamee2 = [0.0; APPROX_LEN];

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
        e1_swamee[i] = e1::swammee_ohija(x).unwrap();
        ei_swamee[i] = ei::convergent_ramanujan(x, 64);
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
    println!("{:#?}", e1_swamee);
    println!("{:#?}", ei_swamee);
    println!("{:#?}", ei_swamee2);


    let node_count = 16_u64;
    let positive = cheb::gen_chebyshev_approximator(10.0,100.0,node_count,(|x| ei::convergent_ramanujan(x, 1024)));
    let negative = cheb::gen_chebyshev_approximator(-10.0,-100.0,node_count,(|x| ei::convergent_ramanujan(x, 1024)));
    let middle_positive = cheb::gen_chebyshev_approximator(0.01,10.0,node_count,(|x| ei::convergent_ramanujan(x, 1024)));
    let middle_negative = cheb::gen_chebyshev_approximator(-0.01,-10.0,node_count,(|x| ei::convergent_ramanujan(x, 1024)));

    let multi_approximator = move |x| -> f64{
        if x < -10.0{
            negative(x)
        }else if x < 0.0{
            middle_negative(x)
        }else if x < 10.0{
            middle_positive(x)
        }else{
            positive(x)
        }
    };
    for i in 0..=2000{
        let x = (i -1000) as f64/10.0;
        let result = multi_approximator(x);

        let result_err = (result-x).abs()/x;
        // println!("{:#?},{:#?},{:#?},{:#?},{:#?},{:#?}", x, ei::convergent_ramanujan(x, 256),pos_err,neg_err,mpos_err,mneg_err);
        println!("{:#?},{:#?},{:#?},{:#?}", x, ei::convergent_ramanujan(x, 256),result,result_err);
    }

    let ret = burmann::integer_partition_array(5);
    println!("{:#?}", ret);

    let ret2 = burmann::integer_partition(5);
    println!("{:#?}", ret2);

    let mut array = vec![0.0; 10];
    array[0] = 1.0;
    for i in 1..10{
        if i & 1 == 0 {
            let n2 = i;
            let n = i/2;
            let sign = if n & 1 == 0{
                1.0
            }else{
                -1.0
            };
            array[i as usize] = 2.0 * sign * (burmann::factorial(n2 - 1) as f64)/(burmann::factorial(n - 1) as f64);
        }else{
            array[i as usize] = 0.0;
        }
    }
    let bc = burmann::calc_burmann_coefficent(1.0, &array, 2, 1);

    let divider = 4.0;
    let node_count = 16;
    let cos_approx = create_cos(node_count, divider);
    let cos_approx2  = create_cos2(node_count, divider);

    const VALUE_LEN:i64 = 17;
    const half_len:i64 = VALUE_LEN/2;
    for i in 0..VALUE_LEN{
        let ix = (i - half_len) as f64;
        let ix = (ix/(divider*half_len as f64));
        let appx = cos_approx(ix * f64::consts::PI);
        let real = f64::cos(ix * f64::consts::PI);
        let appx = cos_approx2(ix * f64::consts::PI);
        println!("#PI*{}, {} vs {} diff {}", ix, appx, real, appx - real);
    }
    // println!("{}", bc);
}


