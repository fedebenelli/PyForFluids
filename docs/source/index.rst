.. PyForFluids documentation master file, created by
   sphinx-quickstart on Mon Nov 15 15:38:39 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyForFluids
===========

|image1| |image2| |image3| |image4|

.. |image1| image:: https://api.codeclimate.com/v1/badges/3551471cd4cdf37e226f/maintainability
   :target: https://codeclimate.com/github/fedebenelli/PyForFluids/maintainability
   :alt: Maintanibility
.. |image2| image:: https://github.com/fedebenelli/pyforfluids/actions/workflows/ci_linux.yml/badge.svg
   :target: https://github.com/fedebenelli/pyforfluids/actions/workflows/ci_linux.yml
   :alt: Tests
.. |image3| image:: https://readthedocs.org/projects/pyforfluids/badge/?version=latest
   :target: https://pyforfluids.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. |image4| image:: https://camo.githubusercontent.com/69644832889fa9dfcdb974614129be2fda8e4591989fd713a983a21e7fd8d1ad/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4469536f6674436f6d7043692d46414d41462d666664613030
   :target: https://camo.githubusercontent.com/69644832889fa9dfcdb974614129be2fda8e4591989fd713a983a21e7fd8d1ad/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4469536f6674436f6d7043692d46414d41462d666664613030
   :alt: DiSoftCompCi

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules
   installation.rst
   tutorial.ipynb

**PyForFluids** (Python-Fortran-Fluids) is a Python package focused in the
calculation of Fluid properties based on Ecuations of State (EoS). It provides
a simple interface to work from Python but also exploits the high performance
Fortran code for the more heavy calculations.

It's designed with modularity in
mind, in a way that new thermodyinamic models are easy to add, they even can be
written either in Python or Fortran.

Available properties
---------------------

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

Motivation
----------

While nowadays there are a lot of tools for calculation of thermodynamic
properties of fluids, most of them either are hard to mantain and don't have an
integrated testing system or are embeded to other softwares (as spredsheat
software) limiting the things that can be done to that enviroment.

PyForFluids aims to be a tool:

* With high performance, since most of it's calculations are done in Fortran
* Easy to scale due to it's modular design using the power of Python objects.
* Continuosly tested (at every `push`)to spot any problems as soon as possible.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
