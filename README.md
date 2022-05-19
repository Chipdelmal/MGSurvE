# MGSurvE: Mosquito Gene SurveillancE

MGSurvE is a project that optimizes mosquito traps' placement in complex heterogeneous landscapes in an effort to minimize the time to detection of genetic variants of interest.

Please have a look at the [documentation](https://chipdelmal.github.io/MGSurvE/) for more info and our [pypi](https://pypi.org/project/MGSurvE/) package for detailed [installation instructions](https://chipdelmal.github.io/MGSurvE/build/html/installation.html), and [tutorials](https://chipdelmal.github.io/MGSurvE/build/html/demos.html).

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/MGSurvE)](https://pypi.org/project/MGSurvE/)
[![PyPI version](https://badge.fury.io/py/MGSurvE.svg)](https://badge.fury.io/py/MGSurvE) 
[![Docker](https://img.shields.io/docker/v/chipdelmal/mgsurve?color=success&label=docker%20img)](https://hub.docker.com/repository/docker/chipdelmal/mgsurve)
[![Unit Tests](https://github.com/chipdelmal/MGSurvE/actions/workflows/PyTests.yml/badge.svg)](https://github.com/Chipdelmal/MGSurvE/tree/main/MGSurvE/test)
[![Flake8](https://github.com/chipdelmal/MGSurvE/actions/workflows/Flake8.yml/badge.svg)](https://github.com/Chipdelmal/MGSurvE/blob/main/.github/workflows/Flake8.yml)
[![Conda](https://github.com/chipdelmal/MGSurvE/actions/workflows/Anaconda.yml/badge.svg)](https://github.com/Chipdelmal/MGSurvE/blob/main/.github/workflows/Anaconda.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Open Source? Yes!](https://badgen.net/badge/Open%20Source%20%3F/Yes%21/blue?icon=github)](https://github.com/Chipdelmal/MGSurvE)


To install the package's latest stable version run:

```
pip install MGSurvE
```

The usage of [anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html) for environment management is strongly recommended.


MGSurvE requires the installation of the [DEAP](https://deap.readthedocs.io/en/master/) optimization package, which should be installed automatically with our previous `pip` command. This package can also be installed with `conda install deap`, if needed; or by having a look at [DEAP's documentation](https://pypi.org/project/deap/) for additional methods,or at our [installation instructions](https://chipdelmal.github.io/MGSurvE/build/html/installation.html) for common issues with some of the dependencies!


Additionally, some of MGSurvE's map-plotting functions use [cartopy](https://scitools.org.uk/cartopy/). Even though the dependency's installation is not strictly required, the easiest way to install this package is with `conda install cartopy`, but in case there are errors in this process, have a look at the package's [installation instructions](https://scitools.org.uk/cartopy/docs/latest/installing.html).



**Under Construction** :construction: **Please check back in a couple of weeks!**

![landscape](https://github.com/Chipdelmal/MGSurvE/raw/main/img/demo.jpg)


# Features and Development

- [x] Different point-types
- [x] Random point-process landscape generators
- [x] Customizable movement kernels
- [x] Customizable trap attractiveness kernels
- [x] Integrated plotting routines
- [x] GA Optimization routines
- [x] Directed network plotting routines
- [x] Demos and Examples
- [x] Immovable traps
- [x] Different attractiveness levels at specific point-type transitions
- [x] Different movement kernels for male and female mosquitos
- [x] Different trap attractiveness levels for male and female mosquitos 
- [x] Thorough timing benchmarks
- [x] Support for xy coordinate systems
- [x] Integrated Map-plotting routines
- [x] Attractiveness/repellency effects

# Authors and Funders

<img src="https://raw.githubusercontent.com/Chipdelmal/pyMSync/master/media/pusheen.jpg" height="125px" align="middle"><img src="https://github.com/Chipdelmal/MGSurvE/blob/main/img/MGSurvE_Logo.png?raw=true" height="125px" align="middle"> <br><br>

* Lead and Dev: [Héctor M. Sánchez C.](https://chipdelmal.github.io/blog/)
* Active Devs: Elijah Bartolome, Lillian Weng, Ayden Salazar, Xingli Yu, Joanna Yoo, Topiltzin Hernandez
* PIs: [David L. Smith](http://www.healthdata.org/about/david-smith), [John M. Marshall](https://publichealth.berkeley.edu/people/john-marshall/)

<br>

<img src="https://github.com/Chipdelmal/MGSurvE/raw/main/img/berkeley.jpg" height="25px"> &nbsp; <img src="https://github.com/Chipdelmal/MGSurvE/raw/main/img/IHME.jpg" height="25px"> &nbsp;  <img src="https://github.com/Chipdelmal/MGSurvE/raw/main/img/IGI.png" height="25px"> &nbsp; <img src="https://github.com/Chipdelmal/MGSurvE/raw/main/img/DARPA.jpg" height="25px"> &nbsp; <img src="https://github.com/Chipdelmal/MGSurvE/raw/main/img/gates.jpg" height="25px">

