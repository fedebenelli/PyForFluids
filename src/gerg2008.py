from . import core
from .fortran import gerg2008f
from .fortran.gerg2008f import thermopropsf

import numpy as np


class GERG2008:

    name = "GERG2008"

    self.valid_components = {
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

    def validate_components(self, components):
        for component in components:
            if component not in self.valid_components:
                raise ValueError(f"{component} ain't  a valid component")

    def validate_ranges(self, temperature, pressure):
        pass

    def calculate_properties(
        self, temperature, pressure, density, composition
    ):
        density_r, temperature_r = gerg2008f.reducing_funcs(x)

        delta = density * density_r
        tau = temperature_r / temperature

        ao = gerg2008f.ideal_term(
            x, density, temperature, density_r, temperature_r
        )
        ar = gerg2008f.residual_term(x, delta, tau)

        z = thermopropsf.zeta(ar)
        cv = thermopropsf.isochoric_heat(tau, r, ao, ar)
        cp = thermopropsf.isobaric_heat(delta, tau, r, ao, ar)
        w = thermopropsf.sound_speed(delta, tau, r, temperature, m, ao, ar)
        isothermal_thermal_coefficent = (
            thermopropsf.isothermal_thermal_coefficent(delta, tau, density, ar)
        )
        dp_dt = thermopropsf.dp_dt(density, delta, tau, r, ar)
        dp_drho = thermopropsf.dp_drho(temperature, delta, r, ar)
        dp_dv = thermopropsf.dp_dv(density, delta, temperature, r, ar)
        p = thermopropsf.pressure(delta, density, r, temperature, ar)
        s = thermopropsf.entropy(tau, r, ao, ar)
        u = thermopropsf.internal_energy(tau, r, temperature, ao, ar)
        h = thermopropsf.enthalpy(delta, tau, r, temperature, ao, ar)
        g = thermopropsf.gibbs_free_energy(delta, r, temperature, ao, ar)
        jt = thermopropsf.joule_thomson_coeff(delta, tau, density, r, ao, ar)
        k = thermopropsf.isentropic_exponent(delta, tau, ao, ar)
        b = thermopropsf.second_thermal_virial_coeff(density_r, ar)
        c = thermopropsf.third_thermal_virial_coeff(density_r, ar)

        return {
            "density_r": density_r,
            "temperature_r": temperature_r,
            "delta": delta,
            "tau": tau,
            "ao": ao,
            "ar": ar,
            "z": z,
            "cv": cv,
            "cp": cp,
            "w": w,
            "isothermal_thermal_coefficent": isothermal_thermal_coefficent,
            "dp_dt": dp_dt,
            "dp_drho": dp_drho,
            "dp_dv": dp_dv,
            "p": p,
            "s": s,
            "u": u,
            "h": h,
            "g": g,
            "jt": jt,
            "k": k,
            "b": b,
            "c": c,
        }
