import pyforfluids as pff
from pyforfluids.equilibrium import flash
import numpy as np
import sys

case = sys.argv[1]

composition = {
    "ethane": 0.02,
    "propane": 0.96,
    "butane": 0.01,
    "isobutane": 0.01
}

temps = {"1": 250, "2": 322.04}

pressure = 1.713e6
gerg = pff.models.GERG2008()

mix = pff.Fluid(
        model=gerg,
        composition=composition,
        temperature=temps[case],
        density=1
)

print("x", "y", "it", "p", "ki")
for i in np.linspace(1, 99, 10):
    z1 = i/100
    z2 = 1 - z1 
    composition = {"ethane": z1, "propane": z2}
    mix.composition = composition
    vapor, liquid, p, it = flash.bub_p(
            mix, temps[case], iterations=100, rtol=1e-5, atol=1e-5
            )
