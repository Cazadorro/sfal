/// Functions that approximate the erfcinv(x), the inverse "Error Function".

use super::erfinv;
use super::consts;
use std::f64;


///https://en.wikipedia.org/wiki/Error_function#Inverse_functions
pub fn maclaurin_series<const MAX_N: usize>(x:f64) -> f64{
    erfinv::maclaurin_series::<MAX_N>(1.0 - x)
}

///https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
pub fn sergei_pade(x:f64)-> f64{
    erfinv::sergei_pade(1.0 - x)
}

///http://www.mimirgames.com/articles/programming/approximations-of-the-inverse-error-function/#
pub fn newton_refinement(x : f64, refinement_count : u64, erf_fn : impl Fn(f64) -> f64) -> f64{
    erfinv::newton_refinement(1.0 - x, refinement_count, erf_fn)
}

///http://www.mimirgames.com/articles/programming/approximations-of-the-inverse-error-function/#
///Numerical Recipes Third Edition p265
pub fn halley_refinement(x : f64, refinement_count : u64, erfc_fn : impl Fn(f64) -> f64) -> f64{
    erfinv::halley_refinement(1.0 - x, refinement_count, erfc_fn)
}


///from "Approximating the erfinv function" by Mike Giles
pub fn mike_giles(x : f64) -> f64{
   erfinv::mike_giles(1.0 - x)
}