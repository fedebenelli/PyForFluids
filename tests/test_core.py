"""
"""
import pyforfluids as pff


def test_fluid_density():
    composition = {"methane": 1}
    temperature = 200
    density = 2

    fluid = pff.Fluid(
        model=pff.models.GERG2008(),
        composition=composition,
        temperature=temperature,
        density=density,
    )

    fluid.set_density(5)

    assert fluid.density == 5


def test_fluid_temperature():
    composition = {"methane": 1}
    temperature = 200
    density = 2

    fluid = pff.Fluid(
        model=pff.models.GERG2008(),
        composition=composition,
        temperature=temperature,
        density=density,
    )

    fluid.set_temperature(5)

    assert fluid.temperature == 5
