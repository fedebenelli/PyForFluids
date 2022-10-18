import numpy as np

import pyforfluids as pff
import pyforfluids.models as models
from pyforfluids.equilibrium.flash import flash_pt


def test_flash_pt(flash):
    composition = flash["z"].to_dict()
    vapor_composition = flash["y"].to_dict()
    liquid_composition = flash["x"].to_dict()

    temperature = 366.48
    pressure = 1.039e6

    fluid = pff.Fluid(
        model=models.GERG2008(),
        composition=composition,
        temperature=temperature,
        density=1,
    )

    vapor, liquid, beta, it = flash_pt(fluid, pressure, temperature)

    print()
    print("feed:", composition)
    print(vapor_composition)
    print(vapor.composition)
    print(liquid_composition)
    print(liquid.composition)

    for component in composition:
        np.testing.assert_allclose(
            vapor_composition[component],
            vapor.composition[component],
            rtol=0.001,
            atol=1e-3,
        )
        np.testing.assert_allclose(
            liquid_composition[component],
            liquid.composition[component],
            rtol=0.001,
            atol=1e-3,
        )
