/******************************************************************************

Special Relativity toolkit

******************************************************************************/
pub mod srtoolkit {
    use crate::constants::*;
    use ndarray::Array1;

    // Lorentz factor
    pub fn v_rela(lorentz_factor: f64) -> f64 {
        if lorentz_factor <= 1.0 {
            0.0
        } else {
            (1.0 - 1.0 / (lorentz_factor.powi(2))).sqrt()
        }
    }

    // gofb_s: Scalar version
    pub fn lorentz_f(beta: f64) -> f64 {
        1.0 / (1.0 - beta.powi(2)).sqrt()
    }

    // pofg_s: Scalar version
    pub fn pofg_s(g: f64) -> f64 {
        (g * g - 1.0).sqrt()
    }

    // pofg_v: Vector version
    pub fn pofg_v(g: &[f64]) -> Vec<f64> {
        g.iter().map(|&gi| (gi * gi - 1.0).sqrt()).collect()
    }

    // gofp_s: Scalar version
    pub fn gofp_s(p: f64) -> f64 {
        (p * p + 1.0).sqrt()
    }

    // gofp_v: Vector version
    pub fn gofp_v(p: &[f64]) -> Vec<f64> {
        p.iter().map(|&pi| (pi * pi + 1.0).sqrt()).collect()
    }


    // ip2ofg_s: Scalar version
    pub fn ip2ofg_s(g: f64) -> f64 {
        let eps = 1e-6;
        if (g - 1.0).powi(3) < 32.0 * eps {
            0.5 / (g - 1.0) - 0.25 + 0.125 * (g - 1.0)
        } else {
            1.0 / (g * g - 1.0)
        }
    }

    // ip2ofg_v: Vector version
    pub fn ip2ofg_v(g: &[f64]) -> Vec<f64> {
        let eps = 1e-6;
        g.iter().map(|&gi| {
            if (gi - 1.0).powi(3) < 32.0 * eps {
                0.5 / (gi - 1.0) - 0.25 + 0.125 * (gi - 1.0)
            } else {
                1.0 / (gi * gi - 1.0)
            }
        }).collect()
    }

    // Doppler
    pub fn doppler(gamma: f64, mu: f64) -> f64 {
        1.0 / (gamma * (1.0 - v_rela(gamma) * mu))
    }

    // mu_obs_f: Scalar version
    pub fn mu_obs_f(gamma: f64, muc: f64) -> f64 {
        (muc + v_rela(gamma)) / (1.0 + v_rela(gamma) * muc)
    }

    // mu_com_f: Scalar version
    pub fn mu_com_f(gamma: f64, muo: f64) -> f64 {
        (muo - v_rela(gamma)) / (1.0 - v_rela(gamma) * muo)
    }

    // nu_obs_s: Scalar version
    pub fn nu_obs_s(nucom: f64, z: f64, dopp: f64) -> f64 {
        nucom * dopp / (1.0 + z)
    }

    // nu_obs_v: Vector version
    pub fn nu_obs_v(nucom: &[f64], z: f64, dopp: f64) -> Vec<f64> {
        nucom.iter().map(|&nuc| nuc * dopp / (1.0 + z)).collect()
    }

    // nu_com_s: Scalar version
    pub fn nu_com_s(nuobs: f64, z: f64, dopp: f64) -> f64 {
        nuobs * (1.0 + z) / dopp
    }

    // nu_com_v: Vector version
    pub fn nu_com_v(nu: &[f64], z: f64, dopp: f64) -> Vec<f64> {
        nu.iter().map(|&n| n * (1.0 + z) / dopp).collect()
    }

    // t_obs_s: Scalar version
    pub fn t_obs_s(tcom: f64, z: f64, gamma: f64, muo: f64) -> f64 {
        (1.0 + z) * tcom / doppler(gamma, muo)
    }

    // t_obs_v: Vector version
    pub fn t_obs_v(tcom: &[f64], z: f64, gamma: f64, muo: f64) -> Vec<f64> {
        tcom.iter().map(|&tc| (1.0 + z) * tc / doppler(gamma, muo)).collect()
    }

    // t_com_s: Scalar version
    pub fn t_com_s(tobs: f64, z: f64, gamma: f64, muo: f64) -> f64 {
        doppler(gamma, muo) * tobs / (1.0 + z)
    }

    // t_com_v: Vector version
    pub fn t_com_v(tobs: &[f64], z: f64, gamma: f64, muo: f64) -> Vec<f64> {
        tobs.iter().map(|&to| doppler(gamma, muo) * to / (1.0 + z)).collect()
    }

    // x_com_s: Scalar version
    pub fn x_com_s(t: f64, tobs: f64, z: f64, gamma: f64, muo: f64) -> f64 {
        CLIGHT * ((t / doppler(gamma, muo)) - (tobs / (1.0 + z))) / (gamma * (muo - v_rela(gamma)))
    }

    // x_com_v: Vector version
    pub fn x_com_v(t: &[f64], tobs: f64, z: f64, gamma: f64, muo: f64) -> Vec<f64> {
        t.iter().map(|&ti| {
            CLIGHT * ((ti / doppler(gamma, muo)) - (tobs / (1.0 + z))) / (gamma * (muo - v_rela(gamma)))
        }).collect()
    }

}
