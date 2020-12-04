#SFAL *Special Function Approximation Library*


##Summary

Crate that implements **multiple approximations** for

- [tgamma] true gamma function (not implimented)
- [lgamma] ln(tgamma(x)) (not implimented)
- [erf] error function
- [erfc] complementary error function
- [erfinv] inverse error function
- [erfcinv] inverse complementary error function
- [ei] the exponential integral
- [e1] the *other* exponential integral

in addition to:

- [cheb] chebyshev polynomials so you have *an additional* automatic method of approximating these functions
 
entirely in rust.  

## Why?

Most crates, or even libraries outside of rust, that implement ei, or e1 don't do it over the whole numeric range (this effectively does).
Most crates either give you the most accurate implementation possible, or a very inaccurate but very fast approximation, or worse, a combination of both. 

This gives you every approximation I could find and at least understand enough to implement here. 

## Future

every single one of these function needs to still be implemented (except expint)
https://en.cppreference.com/w/cpp/numeric/special_functions

need to templatize everything (requires asPrimitive conversions)

need to implement automatic [burmann series generator.](https://www.semanticscholar.org/paper/On-B%C3%BCrmann's-Theorem-and-Its-Application-to-of-and-Sch%C3%B6pf-Supancic/eec2f0f6260e486f8a4fffb7f619e0717fae4645)
currently burmann requires many complicated solving functions pretty much only found in mathematica. 
Burmann series is like a taylor series approximation for f(x), but instead of centering around x-x0, you center around phi(x) - phi(x0) where phi is a analytic function with specific properties (typically the analytic derivative of f(x))