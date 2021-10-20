"""
"""
import numpy as np
import pytest

import pyforfluids as pff
import pyforfluids.models as models


def test_init():
    composition = {"methane": 1}
    temperature = 200
    pressure = 1e6
    density = 2
    model = models.GERG2008()

    with pytest.raises(Exception) as e_info:
        pff.Fluid(
            model=model, composition=composition, temperature=temperature
        )
    print(e_info)

    with pytest.raises(Exception) as e_info:
        pff.Fluid(
            model=model,
            composition=composition,
            temperature=temperature,
            pressure=pressure,
            density=density,
        )
    print(e_info)

    with pytest.raises(Exception) as e_info:
        pff.Fluid(
            model=model,
            composition="gas",
            temperature=temperature,
            density=density,
        )
    print(e_info)


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


def test_properties():
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

    test_props = {
        "density_r": 10.139342719,
        "temperature_r": 190.564,
        "delta": 0.19725144473637551,
        "tau": 0.95282,
        "ao": np.array(
            [
                [-38.29785372, 0.0, 0.0],
                [5.06967136, -51.93969052, 0.0],
                [-25.70156769, -3.36875089, 0.0],
            ]
        ),
        "ar": np.array(
            [
                [-0.20131914, 0.0, 0.0],
                [-0.97974121, -0.42535655, 0.0],
                [0.41371954, -0.29440971, -2.10427854],
            ]
        ),
        "z": 0.8067446316666418,
        "cv": 27.65109691652404,
        "cp": 46.73894882544477,
        "w": 332.1297523714502,
        "isothermal_thermal_coefficent": -0.45478053068788404,
        "dp_dt": 19.991873834260783,
        "dp_drho": 1046.9355622415853,
        "dp_dv": -4187.742248966341,
        "p": 2683062.2604570426,
        "s": -94.74583033772151,
        "u": -82969.22502587165,
        "h": -81627.69389564313,
        "g": -62678.52782809882,
        "jt": 1.073354071716588,
        "k": 1.025331658435987,
        "b": -0.10473585108040255,
        "c": 0.004165029281375034,
    }

    fluid.calculate_properties()

    for prop in test_props:
        test_prop = test_props[prop]
        calc_prop = fluid.properties[prop]
        np.testing.assert_almost_equal(test_prop, calc_prop, 8)
        # assert np.allclose(test_prop, calc_prop)
