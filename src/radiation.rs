use ndarray::{Array, Axis, Array1};
use ndarray::parallel::prelude::*;

use crate::constants::*;
use crate::misc::*;


pub fn syn_emissivity(freq: f64, gg: &Array1<f64>, nn: &Array1<f64>, b: f64,  rma_func: Option<fn(f64, f64) -> f64>) -> f64 {
    let ng = gg.len();
    let mut jnu = 0f64;

    for k in 0..ng-1 {
        if nn[k] > 1e-100 && nn[k + 1] > 1e-100 {
            let mut qq = -((nn[k + 1] / nn[k]).ln()) / ((gg[k + 1] / gg[k]).ln());
            if qq > 8.0 {
                qq = 8.0;
            } else if qq < -8.0 {
                qq = -8.0;
            }
            jnu += j_mb(freq, b, nn[k], gg[k], gg[k + 1], qq, rma_func);
        }
    }

    if jnu < 1e-200 {
        jnu = 0.0;
    }

    jnu
}

pub fn syn_emissivity_full(freqs: &Array1<f64>,g: &Array1<f64>,n: &Array1<f64>,b: f64,with_abs: bool) -> (Array1<f64>, Array1<f64>) {
    let numdf = freqs.len();
    let mut jmbs = Array1::<f64>::zeros(numdf);
    let mut ambs = Array1::<f64>::zeros(numdf);


    let results: Vec<_> = freqs.axis_iter(Axis(0))
        .into_par_iter()
        .map(|freq| {
            let freq = *freq.first().unwrap();
            let jmb = syn_emissivity(freq, &g, &n, b, Some(rma_new));
            let amb = if with_abs { syn_absorption(freq, &g, &n, b, Some(rma_new)) } else { 0.0 };
            (jmb, amb)
        })
        .collect();

    for (i, (jmb, amb)) in results.into_iter().enumerate() {
        jmbs[i] = jmb;
        if with_abs {
            ambs[i] = amb;
        }
    }

    (jmbs, ambs)
}





pub fn j_mb(nu: f64, b: f64, n0: f64, gmin: f64, gmax: f64, qq: f64, rma_func: Option<fn(f64, f64) -> f64>) -> f64 {
    let nu_b = NUCONST * b;
    let chi = nu / nu_b;
    let i2 = rma_qromb(chi, qq, gmin.ln(), gmax.ln(), rma_func);
    let emiss = JMBCONST * nu_b * n0 * i2 * gmin.powf(qq);

    emiss
}

pub fn rma_qromb(chi: f64, q: f64, lga: f64, lgb: f64, rma_func: Option<fn(f64, f64) -> f64>) -> f64{
    const JMAX:usize =60;
    const JMAXP:usize =61;
    const K:usize = 10;
    const KM:usize = 9;
    const EPS:f64 =1e-5;
    let mut h: [f64; JMAX + 1] = [0.0; JMAX + 1];
    let mut s: [f64; JMAX + 1] = [0.0; JMAX + 1];
    let mut qromb = 0.0;
    let mut dqromb = 0.0;

    h[0] = 1.0;

    for j in 0..JMAX  {
        rma_trapzd(chi, q, lga, lgb, &mut s[j], j + 1, rma_func);
        if j >= K - 1 {
            let h_slice = &h[(j - KM)..=j];
            let s_slice = &s[(j - KM)..=j];
            let pol_r = polint(h_slice, s_slice, 0.0);
            (qromb,dqromb) = pol_r.unwrap();
            if dqromb.abs() <= EPS * qromb.abs() {
                return qromb;
            }
        }
        if j < JMAX - 1 {
            s[j + 1] = s[j];
            h[j + 1] = 0.25 * h[j];
        }
    }

    println!("RMA_qromb error");
    println!("chi    = {}", chi);
    println!("q      = {}", q);
    println!("ga     = {}", lga.exp());
    println!("gb     = {}", lgb.exp());
    println!("qromb  = {}", qromb);
    println!("dqromb = {}", dqromb);
    eprint!("RMA_qromb: too many steps");

    qromb
}


pub fn rma_trapzd(chi: f64, q: f64, lga: f64, lgb: f64, s: &mut f64, n: usize, rma_func: Option<fn(f64, f64) -> f64>){
    //testing but this should work to make the rma function optional
    let func = rma_func.unwrap_or_else(|| rma_new);

    if n == 1 {
        let ega = lga.exp();
        let egb = lgb.exp();
        let fa = ega.powf(1.0 - q) * func(chi, ega);
        let fb = egb.powf(1.0 - q) * func(chi, egb);
        *s = 0.5 * (lgb - lga) * (fa + fb);
    } else {
        let it = 2usize.pow(n as u32 - 2);
        let del = (lgb - lga) / it as f64;
        let mut lg = lga + 0.5 * del;
        let mut fsum = 0.0;
        for _i in 0..it {
            let eg = lg.exp();
            fsum += eg.powf(1.0 - q) * func(chi, eg);
            lg += del;
        }
        *s = 0.5 * (*s + del * fsum);
    }
}

pub fn rma_new(chi: f64, g: f64) -> f64 {
    let c1: f64 = 3.2180900500625734e-4;
    let c2: f64 = 6.50532122717873e-1;
    let c3: f64 = 1.5579904689804556e1;
    let mut res: f64;

    if chi > 0.75 / g {
        let x = 2.0 * chi / (3.0 * g.powi(2));
        if x < c1 {
            res = 1.8084180211028020864 * x.powf(1.0 / 3.0);
        } else if x >= c1 && x <= c2 {
            res = (-0.7871626401625178
                - 0.7050933708504841 * x.ln()
                - 0.35531869295610624 * x.ln().powi(2)
                - 0.06503312461868385 * x.ln().powi(3)
                - 0.0060901233982264096 * x.ln().powi(4)
                - 0.00022764616638053332 * x.ln().powi(5)).exp();
        } else if x > c2 && x <= c3 {
            res = (-0.8236455154570651
                - 0.831668613094906 * x.ln()
                - 0.525630345887699 * x.ln().powi(2)
                - 0.22039314697105414 * x.ln().powi(3)
                + 0.01669179529512499 * x.ln().powi(4)
                - 0.028650695862677572 * x.ln().powi(5)).exp();
        } else {
            res = 0.5 * PI * (-x).exp() * (1.0 - 11.0 / (18.0 * x));
        }
    } else {
        res = 0.0;
    }

    res
}


pub fn syn_absorption(freq: f64,  gg: &Array1<f64>, nn: &Array1<f64>, b: f64,rma_func: Option<fn(f64, f64) -> f64>) -> f64 {
    let ng = gg.len();
    let mut anu = 0.0;

    for k in 0..ng-1 {
        if nn[k] > 1e-100 && nn[k + 1] > 1e-100 {
            let mut qq = -(nn[k + 1] / nn[k]).ln() / (gg[k + 1] / gg[k]).ln();
            if qq > 8.0 {
                qq = 8.0;
            } else if qq < -8.0 {
                qq = -8.0;
            }
            anu += a_mb(freq, b, nn[k], gg[k], gg[k + 1], qq, rma_func);
        }
    }

    if anu < 1e-200 {
        anu = 0.0;
    }

    anu
}

pub fn a_mb(nu: f64,b: f64,n0: f64,gmin: f64,gmax: f64,qq: f64, rma_func: Option<fn(f64, f64) -> f64>) -> f64{
    let nu_b = NUCONST * b;
    let chi = nu / nu_b;
    let a2 = arma_qromb(chi, qq, gmin.ln(), gmax.ln(), rma_func);
    let absor = AMBCONST * nu_b * n0 * a2 * gmin.powf(qq) / nu.powi(2);
    absor
}

pub fn arma_qromb(chi: f64, q: f64, lga: f64, lgb: f64, rma_func: Option<fn(f64, f64) -> f64>) -> f64 {
    let (JMAX,JMAXP,K,KM,EPS)=(60,61,10,9,1e-5);
    let mut h: Vec<f64> = vec![0.0; JMAX + 1];
    let mut s: Vec<f64> = vec![0.0; JMAX + 1];
    let mut dqromb = 0.0;
    let mut qromb = 0.0;

    h[0] = 1.0;

    for j in 0..JMAX {
        arma_trapzd(chi, q, lga, lgb, &mut s[j], j + 1, rma_func);
        if j >= K {
            let h_slice = &h[(j - KM)..=j];
            let s_slice = &s[(j - KM)..=j];
            let pol_r = polint(h_slice, s_slice, 0.0);
            (qromb,dqromb) = pol_r.unwrap();
            if dqromb.abs() <= EPS * qromb.abs() {
                return qromb;
            }
        }
        if j < JMAX - 1 {
            s[j + 1] = s[j];
            h[j + 1] = 0.25 * h[j];
        }
    }

    eprintln!("ARMA_qromb error");
    eprintln!("chi    = {}", chi);
    eprintln!("q      = {}", q);
    eprintln!("ga     = {}", lga.exp());
    eprintln!("gb     = {}", lgb.exp());
    eprintln!("qromb  = {}", qromb);
    eprintln!("dqromb = {}", dqromb);
    eprintln!("ARMA_qromb: too many steps");

    qromb
}



pub fn arma_trapzd(chi: f64,q: f64,lga: f64,lgb: f64,s: &mut f64,n: usize, rma_func: Option<fn(f64, f64) -> f64>) {
    let mut del: f64;
    let mut fsum: f64 = 0.0;
    let mut lg: f64;
    let mut eg: f64;

    let func = rma_func.unwrap_or_else(|| rma_new);

    if n == 1 {
        let ega = lga.exp();
        let egb = lgb.exp();
        let fa = ega.powf(-q) * func(chi, ega) * (q + 1.0 + ega.powi(2) / (ega.powi(2) - 1.0));
        let fb = egb.powf(-q) * func(chi, egb) * (q + 1.0 + egb.powi(2) / (egb.powi(2) - 1.0));
        *s = 0.5 * (lgb - lga) * (fa + fb);
    } else {
        let it = 2usize.pow(n as u32 - 2);
        del = (lgb - lga) / it as f64;
        lg = lga + 0.5 * del;
        let mut fsum = 0.0;
        for _i in 0..it {
            eg = lg.exp();
            fsum += eg.powf(-q) * func(chi, eg) * (q + 1.0 + eg.powi(2) / (eg.powi(2) - 1.0));
            lg += del;
        }

        *s = 0.5 * (*s + del * fsum);
    }
}


///////
/////// radiation transfer
///////


pub fn rad_trans_blob(R: f64, jnu: &Array1<f64>, anu: &Array1<f64>) -> Array1<f64> {
    let nf = jnu.len();
    let mut inu = Array1::<f64>::zeros(nf);

    for j in 0..nf {
        inu[j] = 2.0 * R * jnu[j] * opt_depth_blob(R, anu[j],);
    }

    inu
}


pub fn opt_depth_blob(r: f64, absor: f64) -> f64 {
    let tau = f64::max(1e-100, 2.0 * r * absor);
    if tau <= 1e-10 {
        1.0
    } else {
        let mut u = if tau > 100.0 {
            0.5 - 1.0 / tau.powi(2)
        } else if tau >= 0.01 && tau <= 100.0 {
            0.5 * (1.0 - 2.0 * (1.0 - (1.0 + tau) * (-tau).exp()) / tau.powi(2))
        } else {
            (tau / 3.0) - 0.125 * tau.powi(2)
        };
        u = 3.0 * u / tau;
        u
    }
}
