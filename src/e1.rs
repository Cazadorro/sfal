/// Functions that approximate the E1(x), the "Exponential Integral"., also known as the Theis’ well
/// function (W(x)).

use super::consts;
use super::ei;
use super::errors;
use std::f64;

/// Uses standard converging series, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Convergent_series)
/// Convergent series that approximates E1(x), with max iterations before returning.
/// Series may return before max iterations met.
/// # Arguments
/// * `x` - the x in E1(x)
/// * `max_n` - max iterations before returning
///
/// # Examples
/// ```
/// let e1 = convergent_series(23.0, 128);
/// ```
pub fn convergent_series(x:f64, max_n:u64) -> f64{
    //-Ei(-x) = E1(x)
    -ei::convergent_series(-x, max_n)
}

/// Uses Ramanujan's faster converging series, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Convergent_series)
/// Convergent series that approximates E1(x), with max iterations before returning.
/// Series may return before max iterations met.
/// # Arguments
/// * `x` - the x in E1(x)
/// * `max_n` - max iterations before returning
///
/// # Examples
/// ```
/// let e1 = convergent_ramanujan(23.0, 128);
/// ```
pub fn convergent_ramanujan(x:f64, max_n:u32) -> f64{
    //-Ei(-x) = E1(x)
    -ei::convergent_ramanujan(-x, max_n)
}

/// Uses standard diverging series, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Asymptotic_(divergent)_series)
/// Divergent series that approximates E1(x), with max iterations before returning. Works better on larger values.
/// Series may return before max iterations met.
/// # Arguments
/// * `x` - the x in E1(x)
/// * `max_n` - max iterations before returning
///
/// # Examples
/// ```
/// let e1 = divergent_series(23.0, 128);
/// ```
pub fn divergent_series(x:f64, max_n:u64) -> f64{
    //-Ei(-x) = E1(x)
    -ei::divergent_series(-x, max_n)
}


/// Uses standard continued fraction, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Approximations)
/// Continued fraction that approximates E1(x), with max_n iterations. Works better on larger values,
/// will oscillate due to floating point properties, large N doesn't follow a pattern for better result.
/// # Arguments
/// * `x` - the x in E1(x)
/// * `max_n` - number of iterations before returning
///
/// # Examples
/// ```
/// let e1 = continued_fraction(23.0, 19);
/// ```
pub fn continued_fraction(x:f64, max_n:u64) -> f64{
    let mut fraction = x;
    for i in (1..=max_n).rev(){
        fraction = x + (i as f64/(1.0+(i as f64/fraction)));
    }
    f64::exp(-x) / fraction
}

/// Uses Swamee Ohija appoximation, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Approximations)
/// also see "Pump Test Analysis of confined aquifer" by Swamee and Ohija, the apparent source of
/// the approximation I can't get access to.
/// innaccurate compared to others, but does not rely on iteration of any sort.
/// # Arguments
/// * `x` - the x in E1(x), x must be >= 0.0
///
/// # Examples
/// ```
/// let e1 = swammee_ohija(23.0);
/// ```
pub fn swammee_ohija(x:f64) -> Result<f64, errors::DomainError>{
    if x < 0.0{
        Err(errors::DomainError {
            value: x,
            min: 0.0,
            max: std::f64::INFINITY,
        })
    }
    //TODO find a way to make this work for Ei, negative range not valid.
    let A = f64::ln(((0.56146/x + 0.65) * (1.0 + x)));
    let B = x.powi(4)* f64::exp(7.7*x) * (2.0+x).powf(3.7);
    Ok((A.powf(-7.7) + B).powf(-0.13))
}


/// Uses Allen and Hastings approximation, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Approximations)
/// see "Analytical approximation" by Allen, and "Approximation for Digital Computers" by Hastings
/// the approximation I can't get access to.
/// innaccurate compared to others, but does not rely on iteration of any sort.
/// # Arguments
/// * `x` - the x in E1(x)
///
/// # Examples
/// ```
/// let e1 = allen_hastings(23.0);
/// ```
pub fn allen_hastings(x:f64) -> f64{
    let a = [-0.57722, 0.99999, -0.24991, 0.05519, -0.00976, 0.00108];
    let b = [0.26777, 8.63476, 18.05902, 8.57333];
    let c = [3.95850, 21.09965, 25.63296, 9.57332];
    let x3 = [1.0, x, x.powi(2), x.powi(3)];
    let x5 = [1.0, x, x.powi(2), x.powi(3), x.powi(4), x.powi(5)];
    if x <= 1.0{
        let aTx5 = a.iter().zip(&x5).map(|(ai, xi)| (*ai)*(*xi)).sum::<f64>();
        -(x.abs().ln()) + aTx5
    }else{
        let bTx3 = (b.iter().zip(&x3).map(|(bi, xi)|  (*bi)*(*xi)).sum::<f64>());
        let cTx3 = (c.iter().zip(&x3).map(|(ci, xi)| (*ci)*(*xi)).sum::<f64>());
        (f64::exp(-x)/x) * (bTx3/cTx3)
    }
}

/// Uses approximation by Barry et al, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Approximations)
/// see "Approximation for the exponential integral (Theis well function)" by Barry et al.
/// innaccurate compared to others, but does not rely on iteration of any sort.
/// # Arguments
/// * `x` - the x in E1(x), x must be >= 0.0
///
/// # Examples
/// ```
/// let e1 = barry_et_al(23.0);
/// ```
pub fn barry_et_al(x:f64) -> Result<f64, errors::DomainError>{
    //TODO find a way to make this work for Ei, negative range not valid. 100% possible
    //based on interpolating between standard divergent and convergent series for large and small values.
    if x < 0.0{
        Err(errors::DomainError {
            value: x,
            min: 0.0,
            max: std::f64::INFINITY,
        })
    }
    let q = (20.0/47.0) * x.powf((31.0_f64/26.0).sqrt());
    let G = (-consts::EULER_MASCHERONI).exp();
    let b = ((2.0*(1.0 - G))/(G*(2.0 - G))).sqrt();
    let h_inf = ((1.0 - G) * (G*G - 6.0*G + 12.0))/(3.0*G * (2.0 - G).powi(2)*b);
    let h = (1.0/(1.0 + x * x.sqrt())) + (h_inf*q)/(1.0 +q);
    Ok(((-x).exp()/(G + (1.0 - G) * f64::exp(-x/(1.0-G)))) * f64::ln((1.0 + G/x - (1.0 - G)/(h + b*x).powi(2)).abs()))
}