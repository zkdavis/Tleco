# Tleco

`Tleco` stands for both _in the fire_ and _rise_ in the nahuatl language. `Tleco` is a numerical code that simulates particles _in relativistic plasma_, and the _rise of radiation_ from accelerating particles.
### Featured Publications:
Tleco has already been used for several publications (previously called Paramo):

- Rueda-Becerril, J. M., Harrison, A. O., & Giannios, D. 2021, MNRAS, 501, 4092, doi: 10.1093/mnras/staa3925
- Davis, Z., Rueda-Becerril, J. M., & Giannios, D. 2022, MNRAS, 513, 5766, doi: 10.1093/mnras/stac1282
- Combi, L., & Siegel, D. M. 2023, ApJ, 944, 28, doi: 10.3847/1538-4357/acac29

## Installation

Ensure you have Python 3.10 or later installed on your system. Tleco can be installed via pip and works on all major operating systems.
Using and building with Tleco does not require any compiler. A rust compiler is necessary to develop Tleco.

```console
pip install tleco
```
Tleco consist of both Rust functions and Python functions. Follow the example below to get access to correct functions:
```python
import tleco #tleco.rust_functions
from tleco import constants as cons #from tleco gives access to the python files and functions
```
# Examples and Development

Getting started with Tleco is straightforward and several examples can found in [Examples](https://github.com/zkdavis/Tleco/blob/master/examples)].
Follow these instructions to run the example 'broken_pwl_law_full_example.py'.
```console
git clone https://github.com/zkdavis/Tleco.git
cd Tleco
pip install -r requirements.txt
pip install tleco
python examples/broken_pwl_law_full_example.py
```


### Setting up a developer environment

Contributors and developers are welcome to enhance Tleco's features or address issues. Following the steps above 
for running the broken power-law example will get you nearly ready for development. You will still need to install a
rust compiler(https://www.rust-lang.org/tools/install). Once installed if any of the rust files are changed make sure to
rebuild the object with maturin:
```console
maturin develop
```



### Submitting Changes

* Fork the repository.
* Create your feature branch (git checkout -b feature/AmazingFeature).
* Commit your changes (git commit -am 'Add some AmazingFeature').
* Push to the branch (git push origin feature/AmazingFeature).
* Open a pull request.


# Functions
This section is not complete and will be updated over time
### Rust Functions

- **ic_iso_monochrome** - computes the emissivity($\frac{ergs}{cm^3 Sr}$) at frequency nuout(Hz) from inverse Compton (IC) scattering in an isotropic photon field assuming the photonfield is monochromatic [see in [src/radiation.rs:L356](https://github.com/zkdavis/Tleco/blob/master/src/radiation.rs#L356)]
  - **Parameters:**
    - `nuout` (*f64*): frequency(Hz) in the comoving frame to compute emission at.
    - `uext` (*f64*): energy density($\frac{ergs}{cm^-3}$) of the external photon field in the comoving frame.
    - `nuext` (*f64*): frequency(Hz) in the comoving frame of the external photon field.
    - `n` (*&Array1<f64>*): particle distribution as function of lorentz factor
    - `g` (*&Array1<f64>*): Lorentz factor grid
  - **Returns:**
    - `jnu` (*f64*): emissivity($\frac{ergs}{cm^3 Sr}$) for frequency nuout


- **ic_iso_monochrome_full** - computes the emissivity($\frac{ergs}{cm^3 Sr}$) from inverse Compton (IC) scattering in an isotropic photon field assuming the photonfield is monochromatic [see in [src/radiation.rs:L406](https://github.com/zkdavis/Tleco/blob/master/src/radiation.rs#L406)]
  - **Parameters:**
    - `freqs` (*&Array1<f64>*): frequency(Hz) array in the comoving frame to compute emission over.
    - `uext` (*f64*): energy density($\frac{ergs}{cm^-3}$) of the external photon field in the comoving frame.
    - `nuext` (*f64*): frequency(Hz) in the comoving frame of the external photon field.
    - `n` (*&Array1<f64>*): particle distribution as function of lorentz factor
    - `g` (*&Array1<f64>*): Lorentz factor grid
  - **Returns:**
    - `jic` (*Array1<f64>*): emissivity($\frac{ergs}{cm^3 Sr}$) for frequency range freq


- **rad_cool_pwl** - computes the radiative inverse Compton cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$]) from isotropic photon field uu($\frac{ergs}{cm^-3}$) [see in [src/radiation.rs:L502](https://github.com/zkdavis/Tleco/blob/master/src/radiation.rs#L502)]
  - **Parameters:**
    - `gg` (*&Array1<f64>*): Lorentz factor grid
    - `freqs` (*&Array1<f64>*): frequency(Hz) array in the comoving frame.
    - `uu` (*&Array1<f64>*): energy density($\frac{ergs}{cm^-3}$) of photon field in the comoving frame for every frequency in freqs.
    - `with_kn` (*bool*): bool that will include Klein Nishina affects to the cross section when true.
  - **Returns:**
    - `dotg` (*Array1<f64>*): radiative cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$])


- **rad_cool_mono** - computes the radiative inverse Compton cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$]) from isotropic monotonic photon field u0($\frac{ergs}{cm^-3}$) [see in [src/radiation.rs:L557](https://github.com/zkdavis/Tleco/blob/master/src/radiation.rs#L557)]
  - **Parameters:**
    - `gg` (*&Array1<f64>*): Lorentz factor grid
    - `nu0` (*f64*): frequency(Hz) in the comoving frame of the photon field u0.
    - `u0` (*f64*): energy density($\frac{ergs}{cm^-3}$) of photon field in the comoving frame.
    - `with_kn` (*bool*): bool that will include Klein Nishina affects to the cross section when true.
  - **Returns:**
    - `dotg` (*Array1<f64>*): radiative cooling ($\frac{\partial g}{\partial t}$ [$s^{-1}$])

### Python Functions

- **power_law** - Computes a power law function. [see in [pysrc/tleco/misc_func.py:L33](https://github.com/zkdavis/Tleco/blob/master/pysrc/tleco/misc_func.py#L33)]
  - **Parameters:**
  - **Returns:**
    - `f` (**): Array of output values according to the power law.f = np.zeros_like(g)bounds = (g >= g_min) & (g <= g_max)f[bounds] = np.power(g, p)[bounds]if(normalize):f = f/np.trapz(f,g)return f*n0


# Requirements

The following Python packages are used in this project:

- astropy 6.1.4
- ebltable 0.5.2
- matplotlib 3.9.2
- numpy 2.1.2
- pytest 8.3.3
- scipy 1.14.1
- toml 0.10.2
- tqdm 4.66.5
- maturin 1.5
#### Python Version

This project requires Python ^3.10.
