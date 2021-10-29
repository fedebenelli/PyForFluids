"""
"""

import os

from numpy.distutils.core import Extension, setup

rootdir = os.path.normpath(os.path.join(__file__, os.pardir))
fortrandir = os.path.join(rootdir, "pyforfluids", "fortran")

EXTENSIONS = [
    Extension(
        name="pyforfluids.fortran.gerg2008f",
        sources=[
            os.path.join(fortrandir, "parameters.f95"),
            os.path.join(fortrandir, "gerg.f95"),
        ],
    ),
    Extension(
        name="pyforfluids.fortran.thermo_props",
        sources=[
            os.path.join(fortrandir, "parameters.f95"),
            os.path.join(fortrandir, "thermoprops.f95"),
        ],
    ),
]

PACKAGES = ["pyforfluids", "pyforfluids.models", "pyforfluids.fortran"]

setup(
    name="PyForFluids",
    version="0.0.1",
    description="Library for fluid thermodynamics calculations",
    url="https://github.com/fedebenelli/pyforfluids",
    author="Federico Benelli; María Candelaria Arpajou",
    author_email="federico.benelli@mi.unc.edu.ar",
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
)
