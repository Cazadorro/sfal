/// Functions that approximate the Ei(x), the "Exponential Integral".

use super::consts;
use super::e1;
use std::f64;

/// Uses standard converging series, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Convergent_series)
/// Convergent series that approximates Ei(x), with max iterations before returning.
/// Series may return before max iterations met.
/// # Arguments
/// * `x` - the x in Ei(x)
/// * `max_n` - max iterations before returning
///
/// # Examples
/// ```
/// let ei = convergent_series(23.0, 128);
/// ```
pub fn convergent_series(x:f64, max_n:u64) -> f64{
    let mut total_sum = 0.0;
    let mut n_coefficient = 1.0;
    for n in 1..max_n{
        let n_recip = 1.0/n as f64;
        n_coefficient *= x;
        n_coefficient *= n_recip;
        let prev_sum = total_sum;
        total_sum += n_coefficient*n_recip;
        //comparison against itself basically, shouldn't run into pitfalls with float equality.
        if total_sum == prev_sum{
            break;
        }
    }
    consts::EULER_MASCHERONI + x.abs().ln() + (x / 2.0).exp() * total_sum
}

/// Uses Ramanujan's faster converging series, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Convergent_series)
/// Convergent series that approximates Ei(x), with max iterations before returning.
/// Series may return before max iterations met.
/// # Arguments
/// * `x` - the x in Ei(x)
/// * `max_n` - max iterations before returning
///
/// # Examples
/// ```
/// let ei = convergent_ramanujan(23.0, 128);
/// ```
pub fn convergent_ramanujan(x:f64, max_n:u32) -> f64{
    let mut total_sum = 0.0;
    let mut k_coefficient = 0.0;
    let mut n_coefficient = -1.0/(1.0/2.0);
    for n in 1..max_n{
        if n & 1 == 1{
            let k = ((n as u64 - 1)/2);
            k_coefficient += 1.0/(2*k + 1) as f64;
        }
        n_coefficient *= ((-x)/((n * 2) as f64));
        let prev_sum = total_sum;
        total_sum += n_coefficient * k_coefficient;
        //comparison against itself basically, shouldn't run into pitfalls with float equality.
        if total_sum == prev_sum{
            break;
        }
    }
    consts::EULER_MASCHERONI + x.abs().ln() + (x / 2.0).exp() * total_sum
}

/// Uses standard diverging series, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Asymptotic_(divergent)_series)
/// Divergent series that approximates Ei(x), with max iterations before returning. Works better on larger values.
/// Series may return before max iterations met.
/// # Arguments
/// * `x` - the x in Ei(x)
/// * `max_n` - max iterations before returning
///
/// # Examples
/// ```
/// let ei = divergent_series(23.0, 128);
/// ```
pub fn divergent_series(x:f64, max_n:u64) -> f64{
    //starts off iteration zero.
    let mut nfact_xn = 1.0;
    let mut sum = nfact_xn;
    let mut diff = std::f64::INFINITY;
    let mut first_sum = sum;
    for n in 1..max_n{
        nfact_xn *= (n as f64 / x);
        //nfact_xn /= x;
        let last_sum = sum;
        sum += nfact_xn;
        //creates a "U" of precision to approximation, so we need to figure out when we're "out" of the U
        //see https://en.wikipedia.org/wiki/Asymptotic_expansion#:~:text=In%20mathematics%2C%20an%20asymptotic%20expansion,the%20function%20tends%20towards%20a
        let new_diff = (last_sum - sum).abs();
        if new_diff > diff{
            sum = first_sum;
            break;
        }else{
            first_sum = last_sum;
            diff = new_diff;
        }
    }
    (x.exp()/x)*sum
}

/// Uses standard continued fraction, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Approximations)
/// Continued fraction that approximates Ei(x), with max_n iterations. Works better on larger values,
/// will oscillate due to floating point properties, large N doesn't follow a pattern for better result.
/// # Arguments
/// * `x` - the x in Ei(x)
/// * `max_n` - number of iterations before returning
///
/// # Examples
/// ```
/// let ei = continued_fraction(23.0, 19);
/// ```
pub fn continued_fraction(x:f64, max_n:u64) -> f64{
    let mut fraction = 1.0;
    for i in (1..=max_n).rev(){
        fraction = x + (i as f64/(-1.0+(i as f64/fraction)));
    }
    f64::exp(x) / fraction
}


/// Uses Allen and Hastings approximation, see [here](https://en.wikipedia.org/wiki/Exponential_integral#Approximations)
/// see "Analytical approximation" by Allen, and "Approximation for Digital Computers" by Hastings
/// the approximation I can't get access to.
/// innaccurate compared to others, but does not rely on iteration of any sort.
/// # Arguments
/// * `x` - the x in Ei(x)
///
/// # Examples
/// ```
/// let ei = allen_hastings(23.0);
/// ```
pub fn allen_hastings(x:f64) -> f64{
    //Ei(x) = -E1(-x)
    -e1::allen_hastings(-x)
}

