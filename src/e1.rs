
/// Functions that approximate the Ei(x), the "Exponential Integral".

use super::consts;
use super::ei;
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