# Tleco

`Tleco` stands for both _in the fire_ and _rise_ in the nahuatl language. `Tleco` is a numerical code that simulates particles _in relativistic plasma_, and the _rise of radiation_ from accelerating particles.
### Tleco Paper

- Davis, Z., Rueda-Becerril, J. M., & Giannios, D. 2024, ApJ, 976, 182, doi: 10.3847/1538-4357/ad8bc2

### Featured Publications:
Tleco has already been used for several publications (previously called Paramo):

- Rueda-Becerril, J. M., Harrison, A. O., & Giannios, D. 2021, MNRAS, 501, 4092, doi: 10.1093/mnras/staa3925
- Davis, Z., Rueda-Becerril, J. M., & Giannios, D. 2022, MNRAS, 513, 5766, doi: 10.1093/mnras/stac1282
- Combi, L., & Siegel, D. M. 2023, ApJ, 944, 28, doi: 10.3847/1538-4357/acac29
- Desai, D. D., Haggerty, C. C., Shappee, B. J., Tucker, M. A., Davis, Z., Ashall, C., Chomiuk, L., Gootkin, K., Caprioli, D., Bret, A., & Hakobyan, H. 2025, MNRAS, 541, 2197, doi: 10.1093/mnras/staf1117
- Wu, Z.-F., Guevara-Montoya, S., Beniamini, P., Giannios, D., Groselj, D., & Sironi, L. 2026, arXiv:2601.19135, doi: 10.48550/arXiv.2601.19135

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

- **ic_iso_monochrome** - computes the emissivity($\frac{ergs}{cm^3 Sr}$) at frequency nuout(Hz) from inverse Compton (IC) scattering in an isotropic photon field assuming the photon field is monochromatic [see in [src/radiation.rs:L356](https://github.com/zkdavis/Tleco/blob/master/src/radiation.rs#L356)]
  - **Parameters:**
    - `nuout` (*f64*): frequency(Hz) in the comoving frame to compute emission at.
    - `uext` (*f64*): energy density($\frac{ergs}{cm^-3}$) of the external photon field in the comoving frame.
    - `nuext` (*f64*): frequency(Hz) in the comoving frame of the external photon field.
    - `n` (*&Array1<f64>*): particle distribution as function of lorentz factor
    - `g` (*&Array1<f64>*): Lorentz factor grid
  - **Returns:**
    - `jnu` (*f64*): emissivity($\frac{ergs}{cm^3 Sr}$) for frequency nuout


- **ic_iso_monochrome_full** - computes the emissivity($\frac{ergs}{cm^3 Sr}$) from inverse Compton (IC) scattering in an isotropic photon field assuming the photon field is monochromatic [see in [src/radiation.rs:L406](https://github.com/zkdavis/Tleco/blob/master/src/radiation.rs#L406)]
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

- **power_law** - Computes a power law function. [see in [pysrc/tleco/misc_func.py:L34](https://github.com/zkdavis/Tleco/blob/master/pysrc/tleco/misc_func.py#L34)]
  - **Parameters:**
    - `n0` (*Any*): Normalization constant.
    - `g` (*Any*): Array of input values.
    - `p` (*Any*): Power law index.
    - `g_min` (*Any*): Minimum value of Lorentz factor for the power law application.
    - `g_max` (*Any*): Maximum value of Lorentz factor for the power law application.
    - `normalize` (*Any*) (optional): bool if true normalizes the distribution before returning
  - **Returns:**
    - `f` (*Any*): Array of output values according to the power law.


# Requirements

The following Python packages are used in this project:

- astropy 7.2.0
- ebltable 0.5.2
- pytest 8.3.3
- toml 0.10.2
- tqdm 4.66.5
- maturin 1.5
#### Python Version

This project requires Python ^3.10.
