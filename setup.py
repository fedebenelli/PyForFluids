"""
"""

import os
import platform
import setuptools  # noqa

from numpy.distutils.core import Extension, setup  # noqa

ROOT_DIR = os.path.normpath(os.path.join(__file__, os.pardir))
FORTRAN_DIR = os.path.join(ROOT_DIR, "pyforfluids", "fortran")
REQUIREMENTS = ["numpy>=1.21.2"]
ON_WINDOWS = platform.system() == "Windows"
ON_RTD = os.environ.get("READTHEDOCS") == "True"

if ON_WINDOWS:
    extra_link_args = ["-static", "-static-libgfortran", "-static-libgcc"]
else:
    extra_link_args = []

EXTENSIONS = [
    Extension(
        name="pyforfluids.fortran.gerg2008f",
        sources=[
            os.path.join(FORTRAN_DIR, "parameters.f95"),
            os.path.join(FORTRAN_DIR, "gerg.f95"),
        ],
        extra_link_args=extra_link_args,
    ),
    Extension(
        name="pyforfluids.fortran.thermo_props",
        sources=[
            os.path.join(FORTRAN_DIR, "parameters.f95"),
            os.path.join(FORTRAN_DIR, "thermoprops.f95"),
        ],
        extra_link_args=extra_link_args,
    ),
]

PACKAGES = ["pyforfluids", "pyforfluids.models", "pyforfluids.fortran"]

setup(
    name="PyForFluids",
    version="0.1.0-alpha",
    description="Library for fluid thermodynamics calculations",
    url="https://github.com/fedebenelli/pyforfluids",
    author="Federico Benelli; María Candelaria Arpajou",
    author_email="federico.benelli@mi.unc.edu.ar",
    packages=PACKAGES,
    ext_modules=EXTENSIONS if not ON_RTD else [],
    install_requires=REQUIREMENTS,
)
