"""
"""
import numpy as np

import pyforfluids as pff
import pyforfluids.models as models

import pytest


def test_init():
    composition = {"methane": 1}
    temperature = 200
    pressure = 1e6
    density = 2
    model = models.GERG2008()

    with pytest.raises(ValueError):
        pff.Fluid(
            model=model, composition=composition, temperature=temperature
        )

    with pytest.raises(ValueError):
        pff.Fluid(
            model=model,
            composition=composition,
            temperature=temperature,
            pressure=pressure,
            density=density,
        )

    with pytest.raises(ValueError):
        pff.Fluid(
            model=model,
            composition="gas",
            temperature=temperature,
            density=density,
        )

    fluid = pff.Fluid(
        model=model,
        composition=composition,
        temperature=temperature,
        pressure=pressure,
    )

    test_density = 0.6436568662841562
    calc_density = fluid.density

    np.testing.assert_allclose(calc_density, test_density, 8)


def test_density_iterator():
    composition = {"methane": 1}
    temperature = 200
    density = 2
    model = models.GERG2008()

    fluid = pff.Fluid(
        model=model,
        composition=composition,
        temperature=temperature,
        density=density,
    )

    pressure = 2e6
    test_values = (1.395672933295116, 1994710.425259964, 2)
    calc_values = fluid.density_iterator(pressure)

    np.testing.assert_allclose(test_values, calc_values, 8)

    warning_pressure = 1e15
    with pytest.warns(RuntimeWarning):
        fluid.density_iterator(warning_pressure)


def test_fluid_density():
    composition = {"methane": 1}
    temperature = 200
    density = 2
    model = models.GERG2008()

    fluid = pff.Fluid(
        model=model,
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
    model = models.GERG2008()

    fluid = pff.Fluid(
        model=model,
        composition=composition,
        temperature=temperature,
        density=density,
    )

    fluid.set_temperature(5)

    assert fluid.temperature == 5


def test_fluid_pressure():
    composition = {"methane": 1}
    temperature = 200
    density = 2
    model = models.GERG2008()

    fluid = pff.Fluid(
        model=model,
        composition=composition,
        temperature=temperature,
        density=density,
    )

    fluid.set_pressure(5)

    assert fluid.pressure == 5


def test_fluid_composition():
    composition = {"methane": 1}
    temperature = 200
    density = 2
    model = models.GERG2008()

    fluid = pff.Fluid(
        model=model,
        composition=composition,
        temperature=temperature,
        density=density,
    )

    new_composition = {"ethane": 1}

    fluid.set_composition(new_composition)

    assert fluid.composition == new_composition


def test_isotherm(ahmadi323):
    model = models.GERG2008()
    temperature = 323.3
    composition = {
        "methane": 87.9427,
        "ethane": 6.0000,
        "propane": 2.0430,
        "butane": 0.2998,
        "isobutane": 0.1995,
        "carbon_dioxide": 2.0130,
        "nitrogen": 1.5020,
    }

    for compound in composition:
        composition[compound] = composition[compound] / 100

    fluid = pff.Fluid(
        model=model,
        composition=composition,
        temperature=temperature,
        density=1,
    )

    density_range = ahmadi323["density"].values

    values = fluid.isotherm(density_range)

    for prop in ahmadi323.columns:
        np.testing.assert_allclose(values[prop], ahmadi323[prop], rtol=0.05)


def test_getitem():
    composition = {"methane": 1}
    temperature = 200
    density = 2
    model = models.GERG2008()

    fluid = pff.Fluid(
        model=model,
        composition=composition,
        temperature=temperature,
        density=density,
    )

    np.testing.assert_almost_equal(fluid["pressure"], 2683062.2604570426, 8)


def test_repr():
    composition = {"methane": 0.5, "ethane": 0.5}
    temperature = 200
    density = 1
    model = models.GERG2008()

    fluid = pff.Fluid(
        model=model,
        composition=composition,
        temperature=temperature,
        density=density,
    )

    test_string = (
        "Fluid(model=GERG2008, temperature=200, pressure=1273969.7325, "
        "density=1.0000, composition={'methane': 0.5, 'ethane': 0.5})"
    )

    assert fluid.__repr__() == test_string
