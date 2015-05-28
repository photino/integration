//! Numerical integration methods based on quadrature rules.

/// Gauss–Legendre quadrature.
///
/// ## Example
///
/// ```rust
/// use std::f64;
/// use integration::quadrature;
///
/// fn f(x: f64) -> f64 {
///     x.sin()
/// }
///
/// let y = quadrature::gauss_legendre(f, 0_f64, f64::consts::PI, 1e-10);
///
/// assert!((y - 2.0).abs() < 1e-10);
/// ```
pub fn gauss_legendre(f: fn(f64) -> f64, a: f64, b: f64, eps: f64) -> f64 {
    let x1 = -0.745355992499929898803057889577092078_f64;
    let x2 = 0_f64;
    let x3 = 0.745355992499929898803057889577092078_f64;
    let w1 = 0.555555555555555555555555555555555556_f64;
    let w2 = 0.888888888888888888888888888888888889_f64;
    let w3 = 0.555555555555555555555555555555555556_f64;
    let m = 0.5 * (a + b);
    let d = 0.5 * (b - a);
    let y = d * (w1 * f(d * x1 + m) + w2 * f(d * x2 + m) + w3 * f(d * x3 + m));
    if (y - d * (f(a) + f(b))).abs() > eps {
        gauss_legendre(f, a, m, 0.5 * eps) + gauss_legendre(f, m, b, 0.5 * eps)
    } else {
        y
    }
}

