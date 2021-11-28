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

with open("README.md") as fp:
    LONG_DESCRIPTION = fp.read()

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
    name="pyforfluids",
    version="0.0.1-alpha1",
    description="Fluid's thermodynamic properties",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    author="Federico E. Benelli; M. Candelaria Arpajou",
    author_email="federico.benelli@mi.unc.edu.ar",
    url="https://github.com/fedebenelli/pyforfluids",
    licence="MIT",
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
    ext_modules=EXTENSIONS if not ON_RTD else [],
    install_requires=REQUIREMENTS,
)
