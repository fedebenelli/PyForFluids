# PyForFluids 
<a href="https://codeclimate.com/github/fedebenelli/PyForFluids/maintainability">
<img src="https://api.codeclimate.com/v1/badges/3551471cd4cdf37e226f/maintainability"/></a>
<a href="https://github.com/fedebenelli/pyforfluids/actions/workflows/CI.yml">
<img src="https://github.com/fedebenelli/pyforfluids/actions/workflows/CI.yml/badge.svg">
</a> 
<a href='https://pyforfluids.readthedocs.io/en/latest/?badge=latest'>
<img src='https://readthedocs.org/projects/pyforfluids/badge/?version=latest'
alt='Documentation Status'/></a> <a href="https://github.com/leliel12/diseno_sci_sfw">
<img src="https://camo.githubusercontent.com/69644832889fa9dfcdb974614129be2fda8e4591989fd713a983a21e7fd8d1ad/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4469536f6674436f6d7043692d46414d41462d666664613030"></a>
<a href='https://pypi.org/project/pyforfluids/'>
<img src='https://img.shields.io/pypi/v/pyforfluids'>
</a>

PyForFluids (Python-Fortran-Fluids) is a Python package focused in the
calculation of Fluid properties based on Ecuations of State (EoS). It provides
a simple interface to work from Python but also exploits the high performance
Fortran code for the more heavy calculations.

It’s designed with modularity in mind, in a way that new thermodyinamic models
are easy to add, they even can be written either in Python or Fortran.

- Multifluid equations:
	- GERG-2008 [Paper link](https://pubs.acs.org/doi/10.1021/je300655b)

## Available properties
- Reduced Temperature and Density
- Ideal Helmholtz Energy (Ao)
- Residual Helmholtz Energy (Ar)
- Compresibility Factor (Z)
- Isochoric Heat (Cv)
- Isobaric Heat (Cp)
- Speed of sound (w)
- Isothermal throttling coefficent (δ)
- Pressure derivatives:
	- Temperature
	- Density
	- Volume
- Pressure (P)
- Entropy (S)
- Gibbs Free Energy (G)
- Enthalpy (H)
- Joule-Thompson coefficent
- Isoentropic exponent
- Virial Terms:
	- B
	- C

## Motivation
While nowadays there are a lot of tools for calculation of thermodynamic
properties of fluids, most of them either are hard to mantain and don't have an
integrated testing system or are embeded to other softwares (as spredsheat
software) limiting the things that can be done to that enviroment.

PyForFluids aims to be a tool:

- With high performance, since most of it's calculations are done in Fortran
- Easy to scale due to it's modular design using the power of Python objects.
- Continuosly tested (at every `push`)to spot any problems as soon as possible.

## Instalation
For installing _PyForFluids_ you just need to:

```sh
pip install pyforfluids
```

Make sure to check the requirements first!

### Requirements
Be sure to install `numpy`and a fortran compiler previously, since both are
needed for the compilation of `Fortran` code.

#### NumPy
```sh
pip install numpy
```

#### Fortran Compiler

##### Linux
- **Debian-based** (Debian, Ubuntu, Mint,...)

```sh
sudo apt install gfortran
```

- **Arch-based** (Arch, Manjaro, Garuda, ...)

```sh
sudo pacman -S gfortran
```

##### Windows
We recommended using the Windows Subsystem for Linux 
[WSL](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10)

If WSL ain't being used, the native Windows wheels will be download instead,
so no need to worry!

##### MacOS

```sh
brew install gfortran
```

## Authors
Federico E. Benelli (<a href=federico.benelli@mi.unc.edu.ar>federico.benelli@mi.unc.edu.ar</a>); M. Candelaria
Arpajou
