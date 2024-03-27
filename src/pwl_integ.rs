use crate::constants::*;

// powlaw_integ function
pub fn powlaw_integ(x1: f64, x2: f64, y1: f64, y2: f64) -> f64 {
    let mut q = -(y2 / y1).ln() / (x2 / x1).ln();
    if q > 8.0 {
        q = 8.0;
    }
    if q < -8.0 {
        q = -8.0;
    }
    y1 * x1 * pinteg(x2 / x1, q, 1e-6)
}

// Pinteg function
pub fn pinteg(a: f64, s: f64, eps: f64) -> f64 {
    if ln3(a, eps) * (s - 1.0).powi(2) > 6.0 * eps {
        (1.0 - a.powf(1.0 - s)) / (s - 1.0)
    } else {
        ln1(a, eps) - 0.5 * ln2(a, eps) * (s - 1.0)
    }
}

/// Qinteg function
pub fn qinteg(a: f64, s: f64, eps: f64) -> f64 {
    if 0.125 * ln4(a, eps) * (s - 1.0).powi(2) > eps {
        (1.0 - a.powf(1.0 - s) * (1.0 + (s - 1.0) * ln1(a, eps))) / (s - 1.0).powi(2)
    } else {
        0.5 * ln2(a, eps) - ln3(a, eps) * (s - 1.0) / 3.0
    }
}

// Q2integ function
pub fn q2integ(a: f64, s: f64, eps: f64) -> f64 {
    if 0.1 * ln5(a, eps) * (s - 1.0).powi(2) > eps {
        (2.0 - a.powf(1.0 - s)) * (2.0 + (s - 1.0) * (2.0 + (s - 1.0) * ln1(a, eps))) / (s - 1.0).powi(3)
    } else {
        ln3(a, eps) / 3.0 - 0.25 * (s - 1.0) * ln4(a, eps)
    }
}

// sscR function
pub fn sscr(a: f64, b: f64, c: f64, d: f64, alpha: f64, beta: f64) -> f64 {
    let eps = 1e-9;
    c.powf(beta + alpha + 2.0) * a.powf(alpha + 1.0) * pinteg(d / c, -(beta + alpha + 1.0), eps) * pinteg(b / a, -alpha, eps)
}


// sscS function
pub fn sscs(a: f64, b: f64, c: f64, d: f64, alpha: f64, beta: f64) -> f64 {
    let eps = 1e-9;
    let b_ac = b / (a * c);
    let d_c = d / c;
    if ln3(b / (a * d), eps) * (alpha + 1.0).powi(2) > eps * 6.0 {
        c.powf(beta + 1.0) * (b.powf(alpha + 1.0) * pinteg(d_c, -beta, eps) -
            (a * c).powf(alpha + 1.0) * pinteg(d_c, -(alpha + beta + 1.0), eps)) / (alpha + 1.0)
    } else {
        c.powf(alpha + beta + 2.0) * a.powf(alpha + 1.0) * (ln1(b_ac, eps) *
            pinteg(d_c, -beta, eps) - qinteg(d_c, -beta, eps) + 0.5 * (alpha + 1.0) *
            (ln2(b_ac, eps) * pinteg(d_c, -beta, eps) - 2.0 * ln1(b_ac, eps) *
            qinteg(d_c, -beta, eps) + q2integ(d_c, -beta, eps)))
    }
}


// sscG1ISO function
pub fn sscg1iso(a: f64, b: f64, c: f64, d: f64, alpha: f64, beta: f64) -> f64 {
    sscr(a, b, c, d, alpha, beta) - sscr(a, b, c, d, alpha + 1.0, beta)
}

// sscG2ISO function
pub fn sscg2iso(a: f64, b: f64, c: f64, d: f64, alpha: f64, beta: f64) -> f64 {
    sscs(a, b, c, d, alpha, beta) - sscs(a, b, c, d, alpha + 1.0, beta)
}

// get_g1 function
pub fn get_g1(g2: f64, k: f64, q: f64) -> f64 {
    let ttol = 1e-8;
    let eps = 1e-9;
    let mut xn = 100.0;
    let mut x = -1.0;
    while (((xn - x) / x) as f64).abs() > ttol {
        x = xn;
        let f = x * pinteg(x, q - 1.0, eps) - g2 * k * pinteg(x, q, eps);
        let df = pinteg(x, q - 1.0, eps) + x.powf(2.0 - q) - g2 * k * x.powf(-q);
        xn = x - f / df;
    }
    g2 / xn
}

// LN1 function
pub fn ln1(x: f64, eps: f64) -> f64 {
    if 0.5 * (x - 1.0).powi(2) > eps {
        x.ln()
    } else {
        x - 1.0
    }
}

// LN2 function
pub fn ln2(x: f64, eps: f64) -> f64 {
    if (x - 1.0).powi(3) > eps {
        x.ln().powi(2)
    } else {
        (x - 1.0).powi(2)
    }
}

// LN3 function
pub fn ln3(x: f64, eps: f64) -> f64 {
    if 1.5 * (x - 1.0).powi(4) > eps {
        x.ln().powi(3)
    } else {
        (x - 1.0).powi(3)
    }
}

// LN4 function
pub fn ln4(x: f64, eps: f64) -> f64 {
    if 2.0 * (x - 1.0).powi(5) > eps {
        x.ln().powi(4)
    } else {
        (x - 1.0).powi(4)
    }
}

// LN5 function
pub fn ln5(x: f64, eps: f64) -> f64 {
    if 2.5 * (x - 1.0).powi(6) > eps {
        x.ln().powi(5)
    } else {
        (x - 1.0).powi(5)
    }
}
