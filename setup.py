from numpy.distutils.core import Extension, setup

EXTENSIONS = [
    Extension(
        name="pyforfluids.fortran.gerg2008f",
        sources=[
            "pyforfluids/fortran/parameters.f95",
            "pyforfluids/fortran/gerg.f95",
        ],
    ),
    Extension(
        name="pyforfluids.fortran.thermo_props",
        sources=[
            "pyforfluids/fortran/parameters.f95",
            "pyforfluids/fortran/thermoprops.f95",
        ],
    ),
]

if __name__ == "__main__":
    setup(
        name="pyforfluids",
        description="Library for fluid thermodynamics calculations",
        author="Federico Benelli",
        author_email="federico.benelli@mi.unc.edu.ar",
        packages=["pyforfluids"],
        ext_modules=EXTENSIONS,
    )
