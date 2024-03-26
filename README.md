# Paramo

Paramo (PArticles and RAdiation MOnitor) is a toolkit for modeling relativistic outflows. This code solve the Fokker-Planck equation
and calculates emission self consistently. A description of the code and its core algorithms can be found here (paper link.).

### Featured Publications:
Paramo has already been used for several publications:
- Pub1
- Pub2

## Installation

Ensure you have Python 3.9 or later installed on your system. Paramo can be installed via pip and works on all major operating systems.
Using and building with Paramo does not require any compiler. A rust compiler is necessary to develop Paramo.

```console
pip install paramo
```
# Examples and Development

Getting started with Paramo is straightforward and several examples can found in [Examples](https://github.com/zkdavis/PyParamo/blob/master/examples)].
Follow these instructions to run the example 'broken_pwl_law_full_example.py'.
```console
git clone https://github.com/yourusername/paramo.git
cd paramo
pip install -r requirements.txt
pip isntall paramo
python examples/broken_pwl_law_full_example.py
```


### Setting up a developer environment

Contributors and developers are welcome to enhance Paramo's features or address issues. Following the steps above 
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
# Functions
