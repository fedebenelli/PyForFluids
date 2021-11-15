.. PyForFluids documentation master file, created by
   sphinx-quickstart on Mon Nov 15 15:38:39 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyForFluids
===========

|image1| |image2|

.. |image1| image:: https://api.codeclimate.com/v1/badges/3551471cd4cdf37e226f/maintainability
   :target: https://codeclimate.com/github/fedebenelli/PyForFluids/maintainability
.. |image2| image:: https://github.com/fedebenelli/pyforfluids/actions/workflows/ci_linux.yml/badge.svg
   :target: https://github.com/fedebenelli/pyforfluids/actions/workflows/ci_linux.yml


**PyForFluids** (Python-Fortran-Fluids) is a Python package focused in the
calculation of Fluid properties based on Ecuations of State (EoS). It provides
a simple interface to work from Python but also exploits the high performance
Fortran code for the more heavy calculations.

It's designed with modularity in
mind, in a way that new thermodyinamic models are easy to add, they even can be
written either in Python or Fortran.

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
