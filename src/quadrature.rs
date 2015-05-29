//! Numerical integration methods based on quadrature rules.

/// Gauss–Legendre quadrature.
///
/// ## Example
///
/// ```rust
/// use std::f64;
/// use integration::quadrature;
///
/// let f = &|x: f64| x.sin();
///
/// let y = quadrature::gauss_legendre(f, 0.0, f64::consts::PI, 1e-10);
///
/// assert!((y - 2.0).abs() < 1e-10);
/// ```
pub fn gauss_legendre(f: &Fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    let k = 0.6_f64;
    let x1 = k.sqrt();
    let x2 = 0_f64;
    let x3 = -x1;
    let w1 = 5.0 / 9.0;
    let w2 = 8.0 / 9.0;
    let w3 = w1;
    let m = (a + b) / 2.0;
    let d = (b - a) / 2.0;
    let y = d * (w1 * f(x1.mul_add(d, m)) + w2 * f(x2.mul_add(d, m)) + w3 * f(x3.mul_add(d, m)));
    if (y - d * (f(a) + 4.0 * f(m) + f(b)) / 3.0).abs() > eps {
        gauss_legendre(f, a, m, eps / 2.0) + gauss_legendre(f, m, b, eps / 2.0)
    } else {
        y
    }
}

