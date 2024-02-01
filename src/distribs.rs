use ndarray::{Array1, ArrayView1, s, Zip};
use std::f64::consts::E;
use std::option::Option;
use scilib::math::bessel;
use num::complex::Complex64;
use crate::constants::*;
use crate::SRtoolkit::srtoolkit::*;
use crate::misc::*;
use crate::specialf::*;



#[derive(Debug)]
enum CoolingError {
    MissingArguments(&'static str),
    InvalidCoolType,
}



pub fn eq_59_park1995(t: f64, g: ArrayView1<f64>) -> Array1<f64> {
    let size = g.len();
    let mut gf = Array1::<f64>::zeros(size);

    let d_cof: f64 = 1.0;
    let q: f64 = 3.0;
    let a: f64 = 1.0;
    let x0: f64 = ((100.0f64.powi(2)) - 1.0).sqrt();
    let alpha: f64 = (2.0 - q) / 2.0;
    let theta: f64 = 1.0;
    let t_esc: f64 = 1.0;
    let tau: f64 = d_cof * t;
    let x = g.mapv(|g_i| (g_i.powi(2) - 1.0).sqrt());

    let gfa = (1.0 / alpha.abs()) * (1.0 / (2.0 * tau));
    let gfb = x.mapv(|x_i| x_i.powf((1.0 - q + a) / 2.0)) * x0.powf((1.0 - q - a) / 2.0);
    let etta = (q - 1.0 + a) / (2.0 * alpha);
    let order = etta.abs();
    let gfc = x.mapv(|x_i| bessel::i_nu(order, Complex64::new((x_i.powf(alpha) * x0.powf(alpha)) / (2.0 * tau * alpha.powi(2)), 0.0)).re);
    let gfd = x.mapv(|x_i| (-((x_i.powf(2.0 * alpha)) + x0.powf(2.0 * alpha)) / (4.0 * alpha.powi(2) * tau)).exp());
    let gfe = (-theta * tau).exp();

    gf = gfa * gfb * gfc * gfd * gfe * 4.0 * PI;

    for i in 0..size {
        if gf[i].is_nan() {
            gf[i] = 0.0;
        }
        if gf[i] > 1e200 {
            gf[i] = 1e200;
        }
        if gf[i] < 1e-100 {
            gf[i] = 1e-100;
        }
    }

    gf
}

pub fn fp_findif_difu(dt_in: f64, g: &Array1<f64>, nin: &Array1<f64>, gdot_in: &Array1<f64>,
              din: &Array1<f64>, qin: &Array1<f64>, tesc_in: f64, tlc: f64,
              check_params: Option<bool>) -> Array1<f64> {
    let ng = g.len();
    let mut nout = Array1::<f64>::zeros(ng);

    let check_params2 = check_params.unwrap_or(true);

    if check_params2 {
        let mut check_count = 0;

        // Check scalar values for NaN
        if dt_in.is_nan() { check_count += 1; }
        if tesc_in.is_nan() { check_count += 1; }
        if tlc.is_nan() { check_count += 1; }

        // Check vectors for NaN
        check_count += g.iter().filter(|&&val| val.is_nan()).count();
        check_count += nin.iter().filter(|&&val| val.is_nan()).count();
        check_count += gdot_in.iter().filter(|&&val| val.is_nan()).count();
        check_count += din.iter().filter(|&&val| val.is_nan()).count();
        check_count += qin.iter().filter(|&&val| val.is_nan()).count();

        if check_count > 0 {
            panic!("FP_FinDif_difu: at least one input parameter is not correctly defined");
        }
    }

    let dt = dt_in / tlc;
    let tesc = if tesc_in < 1e100 && tesc_in > 1e-100 { tesc_in / tlc } else { tesc_in };
    let gdot = gdot_in.mapv(|x| x * tlc);
    let dd = din.mapv(|x| x * tlc);
    let qq = qin.mapv(|x| x * tlc);

    let dxp2 = &g.slice(s![1..]) - &g.slice(s![..-1]);
    let dxm2 = Array1::from_iter(dxp2.iter().skip(1).chain(std::iter::once(dxp2.last().unwrap())));
    let dx = (&dxp2 + &dxm2) * 0.5;

    let cc_p2 = (&dd.slice(s![1..]) + &dd.slice(s![..-1])) * 0.25;
    let cc_m2 = cc_p2.clone();
    let bb_p2 = (&dd.slice(s![1..]) - &dd.slice(s![..-1])) / &dx + (&gdot.slice(s![1..]) + &gdot.slice(s![..-1])) * 0.5;
    let bb_m2 = bb_p2.clone();

    let ww_p2 = &bb_p2 / &cc_p2;
    let ww_m2 = &bb_m2 / &cc_m2;

    let zz_p2 = ww_p2.map(|&w| if w * 0.5 > 200.0 { 200.0 } else if w * 0.5 < -200.0 { -200.0 } else { w * 0.5 });
    let zz_m2 = ww_m2.map(|&w| if w * 0.5 > 200.0 { 200.0 } else if w * 0.5 < -200.0 { -200.0 } else { w * 0.5 });

    let yy_p2 = zz_p2.map(|&z| if z.abs() < 0.1 { 1.0 - z.powi(2) / 24.0 + 7.0 * z.powi(4) / 5760.0 - 31.0 * z.powi(6) / 967680.0 } else { z.abs() * (-z.abs()).exp() / (1.0 - (-2.0 * z.abs()).exp()) });
    let yy_m2 = zz_m2.map(|&z| if z.abs() < 0.1 { 1.0 - z.powi(2) / 24.0 + 7.0 * z.powi(4) / 5760.0 - 31.0 * z.powi(6) / 967680.0 } else { z.abs() / (z.abs().exp() - (-z.abs()).exp()) });

    let r = nin + dt * &qq;
    let a = -dt * &cc_m2 * &yy_m2 * (-&zz_m2).mapv(f64::exp) / (&dx * &dxm2);
    let b = 1.0 + dt * (&cc_p2 * &yy_p2 * (-&zz_p2).mapv(f64::exp) / &dxp2 + &cc_m2 * &yy_m2 * &zz_m2.mapv(f64::exp) / &dxm2) / &dx + dt / tesc;
    let c = -dt * &cc_p2 * &yy_p2 * &zz_p2.mapv(f64::exp) / (&dx * &dxp2);

    tridag_ser(&a, &b, &c, &r)
}

