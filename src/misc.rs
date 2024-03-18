// ! *************************************************************************
// !  Miscellaneous functions and routines
// !  ====================================
// !
// !  - polint      : subroutine
// !  - an_error    : subroutine
// !  - arth        : function
// !  - locate      : function
// !  - chebev      : function
// !  - char2int    : function
// !  - char2double : function
// !  - dbesseljn   : function
// !!!TODO: Update list
// !
// ! *************************************************************************

use std::fs::File;
use std::io::{self, BufRead};
use ndarray::Array1;

pub fn an_error(string: &str) {
    panic!("ERROR {}", string);
}

pub fn a_warning(string: &str) {
    println!("Warning message: {}", string);
}

pub fn arth(first: f64, increment: f64, n: usize) -> Vec<f64> {
    let mut arth = Vec::new();
    if n > 0 {
        arth.push(first);
        for k in 1..n {
            arth.push(arth[k - 1] + increment);
        }
    }
    arth
}

pub fn locate(xx: &[f64], x: f64, in_bounds: Option<bool>) -> usize {
    let n = xx.len();
    let ascnd = xx[n - 1] >= xx[0];
    let bounds = in_bounds.unwrap_or(false);

    match xx.binary_search_by(|probe| probe.partial_cmp(&x).unwrap()) {
        Ok(index) => index,
        Err(index) => {
            if index == 0 || (index == n && bounds) {
                return index;
            }
            if ascnd {
                if x < xx[index] { index - 1 } else { index }
            } else {
                if x > xx[index] { index - 1 } else { index }
            }
        }
    }
}

pub fn assert_eq(nn: &[i32], string: &str) -> i32 {
    if nn.windows(2).all(|w| w[0] == w[1]) {
        nn[0]
    } else {
        panic!("ERROR: an assert_eq failed with this tag: {}", string);
    }
}

pub fn char2int(c: &str) -> Result<i32, std::num::ParseIntError> {
    c.parse::<i32>()
}

pub fn char2double(c: &str) -> Result<f64, std::num::ParseFloatError> {
    c.parse::<f64>()
}

pub fn zeros1d(n: usize, small: bool) -> Vec<f64> {
    let fill_value = if small { 1e-200 } else { 0.0 };
    vec![fill_value; n]
}



pub fn count_lines(filename: &str) -> io::Result<usize> {
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);
    Ok(reader.lines().filter_map(Result::ok).count())
}

pub fn poly(cof: &[f64], n: usize, x: f64) -> f64 {
    let mut res = cof[n];
    for i in (0..n).rev() {
        res = res * x + cof[i];
    }
    res
}



pub fn polint(xa: &[f64], ya: &[f64], x: f64) -> Result<(f64, f64), &'static str> {
    if xa.len() != ya.len() {
        return Err("polint: xa and ya have different sizes");
    }
    let n = xa.len();
    let mut c = ya.to_vec();
    let mut d = ya.to_vec();
    let mut ns = 0;
    let mut dy = 0.0;
    let mut y = 0.0;
    let mut ho = vec![0.0; n];
    let mut den = vec![0.0; n];

    for (ho_elem, (xa_elem, &x_elem)) in ho.iter_mut().zip(xa.iter().zip(std::iter::repeat(&x))) {
        *ho_elem = xa_elem - x_elem;
    }

    let mut min_diff = f64::MAX;
    for (i, &ho_elem) in ho.iter().enumerate() {
        let diff = ho_elem.abs();
        if diff < min_diff {
            ns = i;
            min_diff = diff;
            y = ya[i];
        }
    }

    for m in 0..n-1 {
        for j in 0..n-m-1 {
            let hoj = ho[j];
            let hoj_m = ho[j+m+1];
            den[j] = hoj - hoj_m;
            if den[j] == 0.0 {
                return Err("polint: calculation failure");
            }
            den[j] = (c[j+1] - d[j]) / den[j];
            d[j] = hoj_m * den[j];
            c[j] = hoj * den[j];
        }
        if 2 * ns < n - m - 1 {
            dy = c[ns];
        } else {
            dy = d[ns - 1];
            ns -= 1;
        }
        y += dy;
    }
    Ok((y, dy))
}

pub fn linint(x1: f64, x2: f64, x: f64, y1: f64, y2: f64) -> Result<f64, &'static str> {
    let dx = (x2 - x1).abs();
    let dx1 = (x1 - x).abs();
    let dx2 = (x2 - x).abs();

    if dx1 + dx2 > dx {
        return Err("linint: x not between x1 and x2");
    }

    let y = y1 + (y2 - y1) * (x - x1) / dx;
    Ok(y)
}


pub fn chebev(x: f64, coef: &[f64], xmin: f64, xmax: f64) -> Result<f64, &'static str> {
    if (x - xmin) * (x - xmax) > 0.0 {
        return Err("x is not in range in chebev");
    }

    let y = (2.0 * x - xmin - xmax) / (xmax - xmin);
    let y2 = 2.0 * y;
    let mut d = 0.0;
    let mut dd = 0.0;

    for &coeff in coef.iter().rev() {
        let temp = d;
        d = y2 * d - dd + coeff;
        dd = temp;
    }

    Ok(y * d - dd + 0.5 * coef[0])
}

pub fn tridag_ser(a: &Array1<f64>, b: &Array1<f64>, c: &Array1<f64>, r: &Array1<f64>) -> Array1<f64> {
    let n = b.len();
    let mut u = vec![0.0; n];
    let mut gam = vec![0.0; n];
    let mut bet = b[0];
    if bet == 0.0 {
        panic!("tridag_ser: error at code stage 1");
    }
    u[0] = r[0] / bet;

    for j in 1..n {
        gam[j] = c[j - 1] / bet;
        bet = b[j] - a[j - 1] * gam[j];
        if bet == 0.0 {
            panic!("tridag_ser: error at code stage 2");
        }
        u[j] = (r[j] - a[j - 1] * u[j - 1]) / bet;
    }

    for j in (0..n - 1).rev() {
        u[j] -= gam[j + 1] * u[j + 1];
    }

    u.into()

}


pub fn trapzd_w2arg(funci: Option<fn(f64, f64) -> f64>, a: f64, b: f64, s: &mut f64, n: u32, p: f64){

    let func = funci.unwrap();

    if n == 1 {
        let values = vec![a, b];
        let func_results: Vec<f64> = values.iter().map(|&v| func(v,p)).collect();
        *s = 0.5 * (b - a) * func_results.iter().sum::<f64>();
    } else {
        let it = 2_u32.pow(n - 2);
        let del = (b - a) / it as f64;
        let values: Vec<f64> = (0..it).map(|i| a + (i as f64 + 0.5) * del).collect();
        let func_results: Vec<f64> = values.iter().map(|&v| func(v,p)).collect();
        let fsum: f64 = func_results.iter().sum();
        *s = 0.5 * (*s + del * fsum);
    }
}


pub fn qromb_w2arg(func: Option<fn(f64, f64) -> f64>, a: f64, b: f64, p: f64) -> Result<f64, &'static str>{
    const JMAX: usize = 20;
    const K: usize = 5;
    const EPS: f64 = 1e-5;

    let mut h = vec![1.0; JMAX + 1];
    let mut s = vec![0.0; JMAX + 1];


    for j in 1..=JMAX {
        trapzd_w2arg(func, a, b, &mut s[j - 1], j as u32, p);

        if j >= K {
            let (qromb, dqromb) = polint(&h[j - K..j], &s[j - K..j], 0.0)?;
            if dqromb.abs() <= EPS * qromb.abs() {
                return Ok(qromb);
            }
        }

        if j < JMAX {
            s[j] = s[j - 1];
            h[j] = 0.25 * h[j - 1];
        }
    }

    Err("qromb: too many steps")
}

pub fn trapzd_arr<F>(x: &[f64], a: f64, b: f64, f: &[f64], s: &mut f64, n: u32, polint: F)
where
    F: Fn(&[f64], &[f64], f64) -> Result<f64, &'static str>,
{
    if n == 1 {
        let fa = polint(x, f, a).unwrap_or(0.0);
        let fb = polint(x, f, b).unwrap_or(0.0);
        *s = 0.5 * (b - a) * (fa + fb);
    } else {
        let it = 2_u32.pow(n - 2) as usize;
        let del = (b - a) / it as f64;
        let mut xx = a + 0.5 * del;
        let mut fsum = 0.0;

        for _ in 0..it {
            let fx = polint(x, f, xx).unwrap_or(0.0);
            fsum += fx;
            xx += del;
        }

        *s = 0.5 * (*s + del * fsum);
    }
}

pub fn qromb_arr<F, G>(x: &[f64], a: f64, b: f64, f: &[f64], trapzd_arr: F, polint: G) -> Result<f64, &'static str>
where
    F: Fn(&[f64], f64, f64, &[f64], &mut f64, u32),
    G: Fn(&[f64], &[f64], f64) -> Result<(f64, f64), &'static str>,
{
    const JMAX: usize = 30;
    const K: usize = 10;
    const EPS: f64 = 1e-5;

    let mut h = vec![1.0; JMAX + 1];
    let mut s = vec![0.0; JMAX + 1];

    for j in 1..=JMAX {
        trapzd_arr(x, a, b, f, &mut s[j - 1], j as u32);

        if j >= K {
            let (res, dqromb) = polint(&h[j - K..j], &s[j - K..j], 0.0)?;
            if dqromb.abs() <= EPS * res.abs() {
                return Ok(res);
            }
        }

        if j < JMAX {
            s[j] = s[j - 1];
            h[j] = 0.25 * h[j - 1];
        }
    }

    Err("qromb_flux: too many steps")
}

pub fn rk2_arr(y: f64, dydx: &[f64], h: f64) -> f64 {
    let ndum = dydx.len();
    y + 0.5 * h * (dydx[ndum - 2] + dydx[ndum - 1])
}

pub fn rk4<F>(y: &[f64], dydx: &[f64], x: f64, h: f64, derivs: F) -> Vec<f64>
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let ndum = y.len();
    assert_eq!(ndum, dydx.len(), "Mismatched sizes in rk4");

    let hh = h * 0.5;
    let h6 = h / 6.0;
    let xh = x + hh;

    let mut yt = y.iter().zip(dydx).map(|(&yi, &dyi)| yi + hh * dyi).collect::<Vec<_>>();
    let mut dyt = derivs(xh, &yt);

    yt = y.iter().zip(&dyt).map(|(&yi, &dyi)| yi + hh * dyi).collect::<Vec<_>>();
    let mut dym = derivs(xh, &yt);

    yt = y.iter().zip(&dym).map(|(&yi, &dyi)| yi + h * dyi).collect::<Vec<_>>();
    dym = dym.iter().zip(&dyt).map(|(&dyi, &dyti)| dyi + dyti).collect::<Vec<_>>();
    dyt = derivs(x + h, &yt);

    y.iter()
        .zip(dydx)
        .zip(dyt.iter())
        .zip(dym.iter())
        .map(|(((yi, dyi), dyti), dymi)| yi + h6 * (dyi + dyti + 2.0 * dymi))
        .collect()
}

pub fn rkck<F>(
    y: &[f64],
    dydx: &[f64],
    x: f64,
    h: f64,
    yout: &mut [f64],
    yerr: &mut [f64],
    derivs: F,
) where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let n = y.len();
    let mut ak2 = vec![0.0; n];
    let mut ak3 = vec![0.0; n];
    let mut ak4 = vec![0.0; n];
    let mut ak5 = vec![0.0; n];
    let mut ak6 = vec![0.0; n];
    let mut ytemp = vec![0.0; n];

    // Cash-Karp parameters
    let b21 = 0.2;
    let b31 = 3.0 / 40.0;
    let b32 = 9.0 / 40.0;
    let b41 = 0.3;
    let b42 = -0.9;
    let b43 = 1.2;
    let b51 = -11.0 / 54.0;
    let b52 = 2.5;
    let b53 = -70.0 / 27.0;
    let b54 = 35.0 / 27.0;
    let b61 = 1631.0 / 55296.0;
    let b62 = 175.0 / 512.0;
    let b63 = 575.0 / 13824.0;
    let b64 = 44275.0 / 110592.0;
    let b65 = 253.0 / 4096.0;
    let c1 = 37.0 / 378.0;
    let c3 = 250.0 / 621.0;
    let c4 = 125.0 / 594.0;
    let c6 = 512.0 / 1771.0;
    let dc1 = c1 - 2825.0 / 27648.0;
    let dc3 = c3 - 18575.0 / 48384.0;
    let dc4 = c4 - 13525.0 / 55296.0;
    let dc5 = -277.00 / 14336.0;
    let dc6 = c6 - 0.25;

    // First step
    for i in 0..n {
        ytemp[i] = y[i] + b21 * h * dydx[i];
    }
    derivs(x + 0.2 * h, &ytemp, &mut ak2);

    // Second step
    for i in 0..n {
        ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
    }
    derivs(x + 0.3 * h, &ytemp, &mut ak3);

    // Third step
    for i in 0..n {
        ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    }
    derivs(x + 0.6 * h, &ytemp, &mut ak4);

    // Fourth step
    for i in 0..n {
        ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
    }
    derivs(x + h, &ytemp, &mut ak5);

    // Fifth step
    for i in 0..n {
        ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
    }
    derivs(x + 0.875 * h, &ytemp, &mut ak6);

    // Accumulate contributions to yout and yerr
    for i in 0..n {
        yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
        yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
    }
}

pub fn check_isnan_s(x: f64) -> i32 {
    if x.is_nan() {
        1
    } else {
        0
    }
}

pub fn check_isnan_v(x: &[f64]) -> i32 {
    for &value in x {
        if value.is_nan() {
            return 1;
        }
    }
    0
}

pub fn iop() -> [f64; 14] {
    [
        0.9999999999999997, 0.2466405579426905, 0.0147898036344585,
        0.000382699355994036, 5.395676869878828e-6, 4.700912200921704e-8,
        2.733894920915608e-10, 1.115830108455192e-12, 3.301093025084127e-15,
        7.209167098020555e-18, 1.166898488777214e-20, 1.378948246502109e-23,
        1.124884061857506e-26, 5.498556929587117e-30
    ]
}

pub fn ioq() -> [f64; 5] {
    [
        0.4463598170691436, 0.001702205745042606, 2.792125684538934e-6,
        2.369902034785866e-9, 8.965900179621208e-13
    ]
}

pub fn iopp() -> [f64; 5] {
    [
        0.1192273748120670, 0.1947452015979746, 0.0762924181600588,
        0.008474903580801549, 0.0002023821945835647
    ]
}

pub fn ioqq() -> [f64; 6] {
    [
        0.2962898424533095, 0.4866115913196384, 0.1938352806477617,
        0.0226167109340046, 0.0006450448095075585, 1.529835782400450e-6
    ]
}

pub fn i1p() -> [f64; 14] {
    [
        0.5, 0.06090824836578078, 0.00240728857454534, 4.622311145544158e-5,
        5.161743818147913e-7, 3.712362374847555e-9, 1.833983433811517e-11,
        6.493125133990706e-14, 1.693074927497696e-16, 3.299609473102338e-19,
        4.813071975603122e-22, 5.164275442089090e-25, 3.846870021788629e-28,
        1.712948291408736e-31
    ]
}

pub fn i1q() -> [f64; 5] {
    [
        0.4665973211630446, 0.001677754477613006, 2.583049634689725e-6,
        2.045930934253556e-9, 7.166133240195285e-13
    ]
}


pub fn i1pp() -> [f64; 5] {
    [
        0.1286515211317124, 0.1930915272916783, 0.06965689298161343,
        0.007345978783504595, 0.0001963602129240502
    ]
}

pub fn i1qq() -> [f64; 6] {
    [
        0.3309385098860755, 0.4878218424097628, 0.1663088501568696,
        0.01473541892809522, 0.0001964131438571051, -0.0001034524660214173
    ]
}




pub fn bessel_i0(x: f64) -> f64 {
    let iopa = iop();
    let ioppa = iopp();
    let ioqa = ioq();
    let ioqqa = ioqq();

    if x.abs() < 15.0 {
        let y = x * x;
        poly(&iopa, 13, y) / poly(&ioqa, 4, 225.0 - y)
    } else {
        let z = 1.0 - 15.0 / x.abs();
        x.abs().exp() * poly(&ioppa, 4, z) / (poly(&ioqqa, 5, z) * (x.abs()).sqrt())
    }
}

pub fn bessel_i1(x: f64) -> f64 {
    let i1pa = i1p();
    let i1ppa = i1pp();
    let i1qa = i1q();
    let i1qqa = i1qq();

    if x.abs() < 15.0 {
        let y = x * x;
        x * poly(&i1pa, 13, y) / poly(&i1qa, 4, 225.0 - y)
    } else {
        let z = 1.0 - 15.0 / x.abs();
        let result = x.abs().exp() * poly(&i1ppa, 4, z) / (poly(&i1qqa, 5, z) * (x.abs()).sqrt());
        if x < 0.0 {
            -result
        } else {
            result
        }
    }
}

pub fn bessel_in(n: f64, x: f64) -> f64 {
    const ACC: f64 = 200.0;

    if n == 0.0 {
        bessel_i0(x)
    } else if n == 1.0 {
        bessel_i1(x)
    } else {
        let tox = 2.0 / x.abs();
        let mut bip = 0.0;
        let mut ans = 0.0;
        let mut bi = 1.0;
        let mut bim;

        for i in (0..=2 * (n * (ACC * n).sqrt()) as i32).rev() {
            bim = bip + i as f64 * tox * bi;
            bip = bi;
            bi = bim;

            if i == n as i32 {
                ans = bip;
            }
        }
        let result = ans * bessel_i0(x) / bi;
        if x < 0.0 { -result } else { result }
    }
}

pub fn bessel_in_v(n: f64, x: &[f64]) -> Vec<f64> {
    x.iter().map(|&xi| bessel_in(n, xi)).collect()
}
