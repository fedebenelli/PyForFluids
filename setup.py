#!/usr/bin/env python
# -*- coding: utf-8 -*-
# License: MIT License (https://tldrlegal.com/license/mit-license)
# Copyright (c) 2021-2022 Federico Benelli and others.
# All rights reserved.

"""This file is for the distribution of pyforfluids."""

# -> IMPORTS ------------------------------------------------------------------

import os

import setuptools  # noqa

from skbuild import setup

# -----------------------------------------------------------------------------


# -> CONSTANTS ----------------------------------------------------------------

PATH = os.path.normpath(os.path.join(__file__, os.pardir))


PACKAGES = [
    "pyforfluids",
    "pyforfluids.models",
    "pyforfluids.fortran",
    "pyforfluids.equilibrium",
]

with open("README.md") as fp:
    LONG_DESCRIPTION = fp.read()

# -----------------------------------------------------------------------------


# -> REQUIREMENTS -------------------------------------------------------------

REQUIREMENTS = ["numpy>=1.21.2", "pandas>=1.3.5", "scipy>=1.7.3"]

# -----------------------------------------------------------------------------


# -> VERSION ------------------------------------------------------------------

INIT_PATH = os.path.join(PATH, "pyforfluids", "__init__.py")

with open(INIT_PATH, "r") as f:
    for line in f:
        if line.startswith("__version__"):
            VERSION = line.split("=", 1)[-1].replace('"', "").strip()
            break

# -----------------------------------------------------------------------------


# -> SETUP --------------------------------------------------------------------

setup(
    name="pyforfluids",
    version=VERSION,
    description="Fluid's thermodynamic properties",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author="Federico E. Benelli; M. Candelaria Arpajou",
    author_email="federico.benelli@mi.unc.edu.ar",
    url="https://github.com/fedebenelli/pyforfluids",
    license="MIT",
    keywords="Thermodynamic, Fluids, Properties, EoS",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Fortran",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
    ],
    packages=PACKAGES,
    install_requires=REQUIREMENTS,
)

# -----------------------------------------------------------------------------
