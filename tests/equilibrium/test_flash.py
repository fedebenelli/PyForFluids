#!/usr/bin/env python
# -*- coding: utf-8 -*-
# License: MIT License (https://tldrlegal.com/license/mit-license)
# Copyright (c) 2021-2022 Federico Benelli and others.
# All rights reserved.
import numpy as np

import pyforfluids as pff
import pyforfluids.models as models
from pyforfluids.equilibrium.flash import flash_pt, solve_rr

import pytest


@pytest.mark.parametrize(
    "k, z, x, y, beta",
    [
        (
            [2.69667746, 1.2014553, 1.52508674, 0.58173779, 0.2935275],
            [0.01, 0.5, 0.15, 0.2, 0.14],
            [0.00892179, 0.49292686, 0.14459213, 0.20614136, 0.14741818],
            [0.0240592, 0.59222959, 0.22051553, 0.11992022, 0.04327129],
            0.07122802734375,
        ),
        (
            [1.34833873, 0.60072765, 0.76254337, 0.2908689, 0.14676375],
            [0.01, 0.5, 0.15, 0.2, 0.14],
            [0.01, 0.5, 0.15, 0.2, 0.14],
            [0.01, 0.5, 0.15, 0.2, 0.14],
            0,
        ),
        (
            [5.39335492, 2.40291061, 3.05017348, 1.16347559, 0.587055],
            [0.01, 0.5, 0.15, 0.2, 0.14],
            [0.01, 0.5, 0.15, 0.2, 0.14],
            [0.01, 0.5, 0.15, 0.2, 0.14],
            1,
        ),
    ],
)
def test_flash_solve_rr(z, k, x, y, beta):
    z = np.array(z)
    k = np.array(k)

    x_calc, y_calc, beta_calc, _ = solve_rr(z, k)

    rtol = 1e-5

    np.testing.assert_allclose(x_calc, x, rtol=rtol)
    np.testing.assert_allclose(y_calc, y, rtol=rtol)
    np.testing.assert_allclose(beta_calc, beta, rtol=rtol)


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
