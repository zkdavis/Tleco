use ndarray::{Array1, ArrayView1, s, Zip};
use std::f64::consts::E;
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

    let d_cof = 1.0;
    let q: f64 = 3.0;
    let a = 1.0;
    let x0 = ((100.0f64.powi(2)) - 1.0).sqrt();
    let alpha = (2.0 - q) / 2.0;
    let theta = 1.0;
    let t_esc = 1.0;
    let tau = d_cof * t;
    let x = g.mapv(|g_i| (g_i.powi(2) - 1.0).sqrt());

    let gfa = (1.0 / alpha.abs()) * (1.0 / (2.0 * tau));
    let gfb = x.mapv(|x_i| x_i.powf((1.0 - q + a) / 2.0)) * x0.powf((1.0 - q - a) / 2.0);
    let etta = (q - 1.0 + a) / (2.0 * alpha);
    let order = etta.abs();
    let gfc = x.mapv(|x_i| bessel_in(order, (x_i.powf(alpha) * x0.powf(alpha)) / (2.0 * tau * alpha.powi(2))));
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

