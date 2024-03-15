
### syn_emissivity(freq: f64, gamma_bins: &Array1<f64>, n_distrib: &Array1<f64>, b_field: f64,  rma_func: Option<fn(f64, f64) :
#### Description of function: computes the emissivity($\frac{ergs}{cm^3 Sr}$) at frequency nuout(Hz) from inverse Compton (IC) scattering in an isotropic photon field assuming the photon
field is monochromatic

#### Parameters:
None

#### Returns:
- None
### ic_iso_monochrome_full(freqs: &Array1<f64>, uext: f64, nuext: f64, n: &Array1<f64>, g: &Array1<f64>) :
#### Description of function: computes the emissivity($\frac{ergs}{cm^3 Sr}$) from inverse Compton (IC) scattering in an isotropic photon field assuming the photon
field is monochromatic

#### Parameters:
None

#### Returns:
- None
### rad_trans_blob(R: f64, jnu: &Array1<f64>, anu: &Array1<f64>) :
#### Description of function: computes the radiative inverse Compton cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$]) from isotropic photon field uu($\frac{ergs}{cm^-3}$)

#### Parameters:
None

#### Returns:
- None
### rad_cool_mono(gg: &Array1<f64>, nu0: f64, u0: f64, with_kn: bool) :
#### Description of function: computes the radiative inverse Compton cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$]) from isotropic monotonic photon field u0($\frac{ergs}{cm^-3}$)

#### Parameters:
None

#### Returns:
- None