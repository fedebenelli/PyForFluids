#!/usr/bin/env python
# -*- coding: utf-8 -*-
# License: MIT License (https://tldrlegal.com/license/mit-license)
# Copyright (c) 2021-2022 Federico Benelli and others.
# All rights reserved.
"""GERG2008 EoS."""
import warnings

import numpy as np

from ..fortran import fgerg2008
from ..fortran.fthermo_props import thermo_props as tp


class GERG2008:
    """GERG2008 equation of state.

    GERG2008 Equation of state described by O. Kunz and W. Wagner [1]_

    The components must be those of the GERG2008 model.
    This class use imported methods from Fortran subroutines for high speed
    calculation of properties.

    Methods
    -------
    validate_components:
        Check if the components belong to the EOS.
    validate_ranges:
        Check in which range of validity are the temperature and pressure.
    normalize:
        Normalize the composition as molar fractions.
    calculate_properties:
        Calculate the properties.

    References
    ----------
    .. [1] O. Kunz and W. Wagner,
       "The GERG-2008 Wide-Range Equation of State for Natural Gases and
       Other Mixtures: An Expansion of GERG-2004", J. Chem. Eng. Data 2012,
       57, 11, 3032–3091. doi:10.1021/je300655b.
       `<https://pubs.acs.org/doi/10.1021/je300655b>`_
    """

    name = "GERG2008"
    r = fgerg2008.parameters.r
    pc = fgerg2008.parameters.p_c
    tc = fgerg2008.parameters.t_c
    w = fgerg2008.parameters.acentric_factor

    valid_components = {
        "methane",
        "nitrogen",
        "carbon_dioxide",
        "ethane",
        "propane",
        "butane",
        "isobutane",
        "pentane",
        "isopentane",
        "hexane",
        "heptane",
        "octane",
        "nonane",
        "decane",
        "hydrogen",
        "oxygen",
        "carbon_monoxide",
        "water",
        "hydrogen_sulfide",
        "helium",
        "argon",
    }

    names = [
        "methane",
        "nitrogen",
        "carbon_dioxide",
        "ethane",
        "propane",
        "butane",
        "isobutane",
        "pentane",
        "isopentane",
        "hexane",
        "heptane",
        "octane",
        "nonane",
        "decane",
        "hydrogen",
        "oxygen",
        "carbon_monoxide",
        "water",
        "hydrogen_sulfide",
        "helium",
        "argon",
    ]

    def validate_components(self, components):
        """Validate fluid components.

        Verify if the given fluid components are valid
        for the GERG2008 equation of state.If not, a
        ValueError Exception is raised.

        Parameters
        ----------
        components: set
            Set of the fluid components to be verified.
        """
        given_components = set(components)

        diff = given_components.difference(self.valid_components)
        if len(diff) > 0:
            warnings.warn(
                f"{self.name} Valid Components:\n{self.valid_components}",
                category=UserWarning,
            )
            raise ValueError(f"'{diff}' aren't valid components")

    def validate_ranges(self, temperature, pressure):
        """Validate fluid temperature and pressure.

        Verify whether the fluid temperature and pressure values belong to the
        normal, extended or invalid use range of the GERG2008 equation of
        state.  A warning menssage is sent if the temperature and pressure
        conditions are those of the extended or invalid range, and also if they
        take negative values.

        Parameters
        ----------
        temperature: float
            Fluid temperature in Kelvin degrees [K]
        pressure: float
            Fluid pressure in Pascal [Pa]
        """
        if temperature <= 0.0 or pressure <= 0.0:
            warnings.warn(
                "Pressure and Temperature cannot take negative values.",
                category=UserWarning,
            )
        elif temperature > 700.0 or temperature < 60.0:
            warnings.warn(
                "Working conditions belong to the invalid vality range.",
                category=UserWarning,
            )
        elif pressure > 70.0e6:
            warnings.warn(
                "Working conditions belong to the invalid vality range.",
                category=UserWarning,
            )
        elif (temperature >= 90.0 and temperature <= 450.0) and (
            pressure <= 35.0e6
        ):
            return
        else:
            warnings.warn(
                "Working conditions belong to the extended vality range.",
                category=UserWarning,
            )
            return

    def set_concentration(self, composition):
        """
        Verify if the sum of the molar fractions of the fluid components is 1.

        If not, a warninig message is sent and the composition is normalized.

        Parameters
        ----------
        components: dict
            Dictionary of the fluid compounds as keys and their molar
            fraction as values.

        Returns
        -------
        array
            Array of the normalized fluid composition vector.
        """
        methane = composition.get("methane", 0)
        nitrogen = composition.get("nitrogen", 0)
        carbon_dioxide = composition.get("carbon_dioxide", 0)
        ethane = composition.get("ethane", 0)
        propane = composition.get("propane", 0)
        butane = composition.get("butane", 0)
        isobutane = composition.get("isobutane", 0)
        pentane = composition.get("pentane", 0)
        isopentane = composition.get("isopentane", 0)
        hexane = composition.get("hexane", 0)
        heptane = composition.get("heptane", 0)
        octane = composition.get("octane", 0)
        nonane = composition.get("nonane", 0)
        decane = composition.get("decane", 0)
        hydrogen = composition.get("hydrogen", 0)
        oxygen = composition.get("oxygen", 0)
        carbon_monoxide = composition.get("carbon_monoxide", 0)
        water = composition.get("water", 0)
        hydrogen_sulfide = composition.get("hydrogen_sulfide", 0)
        helium = composition.get("helium", 0)
        argon = composition.get("argon", 0)

        x = np.array(
            [
                methane,
                nitrogen,
                carbon_dioxide,
                ethane,
                propane,
                butane,
                isobutane,
                pentane,
                isopentane,
                hexane,
                heptane,
                octane,
                nonane,
                decane,
                hydrogen,
                oxygen,
                carbon_monoxide,
                water,
                hydrogen_sulfide,
                helium,
                argon,
            ],
            dtype="d",
        )

        sum_value = x.sum()

        if sum_value > 1.0001 or sum_value < 0.9999:
            # warnings.warn(
            #    "Composition doesn't add '1', will be normalized",
            #    category=UserWarning,
            # )
            x = x / sum_value

        return x

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
            Dictionary of the thermodynamic properties of the given fluid
            calculated with GERG2008 equation of state.
        """
        # General use parameteres
        fgerg2008.get_params()
        molecular_weights = fgerg2008.parameters.m
        r = self.r

        # Concentration dependant parameters
        x = self.set_concentration(composition)
        m = tp.mean_molecular_weight(x, molecular_weights)
        n = len(x)

        # Reducing functions
        (
            reducing_density,
            reducing_temperature,
            dvr_dx,
            dtr_dx,
            dvr2_dx2,
            dtr2_dx2,
            dvr2_dxx,
            dtr2_dxx,
        ) = fgerg2008.reducing_funcs(x)

        delta = density / reducing_density
        tau = reducing_temperature / temperature

        # Model terms
        ao = fgerg2008.ideal_term(
            x, density, temperature, reducing_density, reducing_temperature
        )

        (
            ar,
            ar_x,
            ar_dx,
            ar_tx,
            ar_xx,
        ) = fgerg2008.residual_term(x, delta, tau)

        if ideal:
            ar = np.zeros_like(ar)
            ar_x = np.zeros_like(ar_x)
            ar_dx = np.zeros_like(ar_dx)
            ar_tx = np.zeros_like(ar_tx)
            ar_xx = np.zeros_like(ar_xx)

        # Properties
        z = tp.zeta(delta, ar)
        p = tp.pressure(delta, density, r, temperature, ar)
        s = tp.entropy(tau, r, ao, ar)
        u = tp.internal_energy(tau, r, temperature, ao, ar)
        h = tp.enthalpy(delta, tau, r, temperature, ao, ar)
        g = tp.gibbs_free_energy(delta, r, temperature, ao, ar)
        jt = tp.joule_thomson_coeff(delta, tau, density, r, ao, ar)
        k = tp.isentropic_exponent(delta, tau, ao, ar)
        cv = tp.isochoric_heat(tau, r, ao, ar)
        cp = tp.isobaric_heat(delta, tau, r, ao, ar)
        w = tp.sound_speed(delta, tau, r, temperature, m, ao, ar)
        isothermal_thermal_coefficent = tp.isothermal_thermal_coefficent(
            delta, tau, density, ar
        )

        # PVT derivatives
        dp_dt = tp.dp_dt(density, delta, tau, r, ar)
        dp_drho = tp.dp_drho(temperature, delta, r, ar)
        dp_dv = tp.dp_dv(density, delta, temperature, r, ar)

        # Equilibrium properties
        (dar_dn, dar_ddn, dar_dtn, dp_dn, dar2_dnn) = tp.molar_derivatives(
            x,
            delta,
            tau,
            r,
            reducing_density,
            reducing_temperature,
            ar,
            ar_x,
            ar_xx,
            ar_tx,
            ar_dx,
            dvr_dx,
            dtr_dx,
            dvr2_dx2,
            dtr2_dx2,
            dvr2_dxx,
            dtr2_dxx,
        )

        msk = np.where(x != 0, 1, 0)
        dnar_dn = (ar[0, 0] + dar_dn) * msk
        dnar2_dtn = -tau / temperature * (ar[1, 1] + dar_dtn)

        dnar2_dnn = dar_dn + dar2_dnn

        excess_volume = -dp_dn / dp_dv

        dp_dn = dp_dn.reshape((1, n))
        dp_dnn = dp_dn.T @ dp_dn

        dlnfug_dt = (
            dnar2_dtn
            + 1 / temperature
            - excess_volume / (r * temperature) * dp_dt
        )
        dlnfug_dp = excess_volume / (r * temperature) - 1 / p
        dlnfug_dn = dnar2_dnn + 1 + dp_dnn / dp_dv / (r * temperature)

        lnfug = dnar_dn - np.log(z)

        # Virial coeficents
        ar_virial, *_ = fgerg2008.residual_term(x, 1e-15, tau)
        b = tp.second_thermal_virial_coeff(reducing_density, ar_virial)
        c = tp.third_thermal_virial_coeff(reducing_density, ar_virial)

        return {
            "critical_density": reducing_density,
            "critical_temperature": reducing_temperature,
            "ideal_helmholtz": ao,
            "residual_helmholtz": ar,
            "ar_x": ar_x,
            "ar_dx": ar_dx,
            "ar_tx": ar_tx,
            "ar_xx": ar_xx,
            "dvr_dx": dvr_dx,
            "dtr_dx": dtr_dx,
            "dvr2_dx2": dvr2_dx2,
            "dtr2_dx2": dtr2_dx2,
            "dvr2_dxx": dvr2_dxx,
            "dtr2_dxx": dtr2_dxx,
            "dar_dn": dar_dn,
            "dadr_dn": dar_ddn,
            "dnar_dn": dnar_dn,
            "dnar2_dtn": dnar2_dtn,
            "dar_ddn": dar_ddn,
            "dar_dtn": dar_dtn,
            "dar2_dnn": dar2_dnn,
            "dnar2_dnn": dar2_dnn,
            "dp_dn": dp_dn,
            "dp_dnn": dp_dnn,
            "excess_volume": excess_volume,
            "compressibility_factor": z,
            "isochoric_heat": cv,
            "isobaric_heat": cp,
            "sound_speed": w,
            "isothermal_thermal_coefficent": isothermal_thermal_coefficent,
            "dp_dt": dp_dt,
            "dp_drho": dp_drho,
            "dp_dv": dp_dv,
            "pressure": p,
            "entropy": s,
            "internal_energy": u,
            "enthalpy": h,
            "gibbs_free_energy": g,
            "joule_thomson_coefficent": jt,
            "isentropic_exponent": k,
            "second_thermal_virial_coeff": b,
            "third_thermal_virial_coeff": c,
            "lnfug": lnfug,
            "dlnfug_dt": dlnfug_dt,
            "dlnfug_dp": dlnfug_dp,
            "dlnfug_dn": dlnfug_dn,
        }

    def __repr__(self):
        """Model Name."""
        return self.name
