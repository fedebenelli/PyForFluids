#!/usr/bin/env python
# -*- coding: utf-8 -*-
# License: MIT License (https://tldrlegal.com/license/mit-license)
# Copyright (c) 2021-2022 Federico Benelli and others.
# All rights reserved.
"""Cubic EoS models."""

import numpy as np

from pyforfluids.fortran.pyfeos import ar


class CubicEOS:
    """Cubic EoS object.

    Represents a cubic equation of state model.

    Attributes
    ----------
    model: str
        Equation of State to be used. Available models:
            - 'PR': PengRobinson Equation of State
            - 'SRK': SRK Equation of State
    mix_rule: str
        Mixing rule of fluids.
            - 'ClassicVdW': classic quadratic mixing rules.
    critical_temperature: array_like
        Array of critical temperatures.
    critical_pressure: array_like
        Array of critical pressures.
    accentric_factor: array_like
        Array of accentric factors.
    kij_matrix: 2D-array_like
        Constant binary interaction parameters.
    lij_matrix: 2D-array_like
        Binary repulsion parameters.

    Methods
    -------
    validate_components:
        For Fluid compatibility, always returns True.
    set_concentration:
        Get a concentrations vector from the composition dictionary.
    calculate_properties:
        Calculate the available properties.
    """

    # Gas constant from inside the Fortran code
    r = 0.08314472

    name = "Cubic EoS"

    __models = {"SRK", "PR"}
    __mixing_rules = {"ClassicVdW"}

    def __init__(
        self,
        model,
        mix_rule,
        names,
        critical_temperature,
        critical_pressure,
        acentric_factor,
        kij_matrix,
        lij_matrix,
    ):
        self.model = model
        self.names = names
        self.tc = np.array(critical_temperature)
        self.pc = np.array(critical_pressure)
        self.w = np.array(acentric_factor)
        self.mix_rule = mix_rule
        self.kij = kij_matrix
        self.lij = lij_matrix

    def validate_components(self, composition):
        """Check if the components are included in the model.

        Parameters
        ----------
        composition: dictionary
            Composition dictionary.

        Returns
        -------
        bool:
            Boolean that represents if the composition is included.
        """
        return True

    def set_concentration(self, composition):
        """Return a composition vector from a composition dictionary.

        Parameters
        ----------
        composition: dict
            Composition dictionary

        Returns
        -------
        array_like:
            Compositions vector
        """
        return np.array([composition[i] for i in composition], dtype="d")

    def calculate_properties(
        self, temperature, pressure, density, composition, ideal=False
    ):
        """Calculate the thermodynamic properties of the given fluid.

        Calculation of the thermodynamic properties of the fluid at it's given
        temperature and density.

        Parameters
        ----------
        temperature: float
            Fluid temperature in Kelvin degrees [K]
        pressure: float
            Fluid pressure in Pascal [Pa]
        density: float
            Fluid density in mol per liter [mol/L]
        composition: dict
            Dictionary of the compounds concentrations as:
            ``{"methane": 0.8, "ethane":0.2}``
            When necessary, the concentration values are normalized.

        Returns
        -------
        dict
            Dictionary of the thermodynamic properties of the given fluid.
        """
        volume = 1 / density
        rt = self.r * temperature

        concentrations = self.set_concentration(composition)
        number_of_moles = concentrations.sum()

        (
            ar_val,
            ar_dt,
            ar_dv,
            ar_dt2,
            ar_dv2,
            ar_dtv,
            ar_dn,
            ar_dn2,
            ar_dtn,
            ar_dvn,
        ) = ar(
            volume,
            temperature,
            self.model,
            self.mix_rule,
            concentrations,
            self.pc,
            self.tc,
            self.w,
            self.kij,
            self.lij,
        )

        # Properties
        p = rt / volume - ar_dv
        z = (p * volume) / (rt)
        dp_dv = -rt * ar_dv2 - rt / volume
        dp_dt = -rt * ar_dtv + p / temperature
        dp_dn = -rt * ar_dvn + rt / volume
        dv_dn = -dp_dn / dp_dv
        dp_drho = ar_dv2 * volume**2 + number_of_moles * rt

        lnfug = ar_dn / (rt) - np.log(z)
        dlnfug_dt = ar_dtn + 1 / temperature - dv_dn / rt * dp_dt
        dlnfug_dp = dv_dn / rt - 1 / p
        dlnfug_dn = (
            number_of_moles * ar_dn2
            + 1
            + number_of_moles
            / rt
            * np.array([dp_dn])
            * np.array([dp_dn]).T
            / dp_dv
        )

        return {
            "compressibility_factor": z,
            "pressure": p * 1e5,
            "dp_dv": dp_dv,
            "dp_dt": dp_dt,
            "dp_dn": dp_dn,
            "dp_drho": dp_drho * 1e2,
            "residual_helmholtz": ar_val,
            "dar_dt": ar_dt,
            "dar_dv": ar_dv,
            "dar_dt2": ar_dt2,
            "dar_dv2": ar_dv2,
            "dar_dtv": ar_dtv,
            "dar_dn": ar_dn,
            "dar_dn2": ar_dn2,
            "dar_dtn": ar_dtn,
            "dar_dvn": ar_dvn,
            "lnfug": lnfug,
            "dlnfug_dt": dlnfug_dt,
            "dlnfug_dp": dlnfug_dp,
            "dlnfug_dn": dlnfug_dn,
        }

    def __check_model(self, model):
        """Check if the model is available.

        Parameters
        ----------
        model: str
            Model to be used, can be checked with:
                `pyforfluids.models.CubicEOS.model_selector`

        Returns
        -------
        bool
            Is the selected model available?
        """
        return model in self.__models

    def __repr__(self):
        """Model representation."""
        return f"{self.name}: {self.model}"
