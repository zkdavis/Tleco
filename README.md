# Functions
### Rust Functions

- **ic_iso_monochrome** - computes the emissivity($\frac{ergs}{cm^3 Sr}$) at frequency nuout(Hz) from inverse Compton (IC) scattering in an isotropic photon field assuming the photonfield is monochromatic [see in [src/radiation.rs:L356](https://github.com/zkdavis/PyParamo/blob/master/src/radiation.rs#L356)]
  - **Parameters:**
    - `nuout` (*f64*): frequency(Hz) in the comoving frame to compute emission at.
    - `uext` (*f64*): energy density($\frac{ergs}{cm^-3}$) of the external photon field in the comoving frame.
    - `nuext` (*f64*): frequency(Hz) in the comoving frame of the external photon field.
    - `n` (*&Array1<f64>*): particle distribution as function of lorentz factor
    - `g` (*&Array1<f64>*): Lorentz factor grid
  - **Returns:**
    - `jnu` (*f64*): emissivity($\frac{ergs}{cm^3 Sr}$) for frequency nuout


- **ic_iso_monochrome_full** - computes the emissivity($\frac{ergs}{cm^3 Sr}$) from inverse Compton (IC) scattering in an isotropic photon field assuming the photonfield is monochromatic [see in [src/radiation.rs:L406](https://github.com/zkdavis/PyParamo/blob/master/src/radiation.rs#L406)]
  - **Parameters:**
    - `freqs` (*&Array1<f64>*): frequency(Hz) array in the comoving frame to compute emission over.
    - `uext` (*f64*): energy density($\frac{ergs}{cm^-3}$) of the external photon field in the comoving frame.
    - `nuext` (*f64*): frequency(Hz) in the comoving frame of the external photon field.
    - `n` (*&Array1<f64>*): particle distribution as function of lorentz factor
    - `g` (*&Array1<f64>*): Lorentz factor grid
  - **Returns:**
    - `jic` (*Array1<f64>*): emissivity($\frac{ergs}{cm^3 Sr}$) for frequency range freq


- **rad_cool_pwl** - computes the radiative inverse Compton cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$]) from isotropic photon field uu($\frac{ergs}{cm^-3}$) [see in [src/radiation.rs:L502](https://github.com/zkdavis/PyParamo/blob/master/src/radiation.rs#L502)]
  - **Parameters:**
    - `gg` (*&Array1<f64>*): Lorentz factor grid
    - `freqs` (*&Array1<f64>*): frequency(Hz) array in the comoving frame.
    - `uu` (*&Array1<f64>*): energy density($\frac{ergs}{cm^-3}$) of photon field in the comoving frame for every frequency in freqs.
    - `with_kn` (*bool*): bool that will include Klein Nishina affects to the cross section when true.
  - **Returns:**
    - `dotg` (*Array1<f64>*): radiative cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$])


- **rad_cool_mono** - computes the radiative inverse Compton cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$]) from isotropic monotonic photon field u0($\frac{ergs}{cm^-3}$) [see in [src/radiation.rs:L557](https://github.com/zkdavis/PyParamo/blob/master/src/radiation.rs#L557)]
  - **Parameters:**
    - `gg` (*&Array1<f64>*): Lorentz factor grid
    - `nu0` (*f64*): frequency(Hz) in the comoving frame of the photon field u0.
    - `u0` (*f64*): energy density($\frac{ergs}{cm^-3}$) of photon field in the comoving frame.
    - `with_kn` (*bool*): bool that will include Klein Nishina affects to the cross section when true.
  - **Returns:**
    - `dotg` (*Array1<f64>*): radiative cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$])

### Python Functions

- **power_law** - Computes a power law function. [see in [pysrc/paramo/misc_func.py:L33](https://github.com/zkdavis/PyParamo/blob/master/pysrc/paramo/misc_func.py#L33)]
  - **Parameters:**
  - **Returns:**
    - `f` (**): Array of output values according to the power law.f = np.zeros_like(g)bounds = (g >= g_min) & (g <= g_max)f[bounds] = np.power(g, p)[bounds]if(normalize):f = f/np.trapz(f,g)return f*n0


# Requirements

The following Python packages are used in this project:

- astropy 6.0.0
- ebltable 0.5.2
- matplotlib 3.8.3
- numpy 1.26.4
- pytest 8.1.1
- toml 0.10.2
- tqdm 4.66.2
- maturin 1.5
#### Python Version

This project requires Python ^3.9.
