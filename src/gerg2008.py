""""""
import warnings

import numpy as np

from fortran import gerg2008f
from fortran.gerg2008f import thermo_props


class GERG2008:
    def __init__(self):
        self.name = "GERG2008"

        self.valid_components = [
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
        for component in components:
            if component not in self.valid_components:
                raise ValueError(f"{component} ain't  a valid component")

    def validate_ranges(self, temperature, pressure):
        pass

    def set_concentration(self, composition):
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
            warnings.warn("Composition doesn't add '1', will be normalized")
            x = x / sum_value

        return x

    def calculate_properties(
        self, temperature, pressure, density, composition
    ):

        gerg2008f.get_params()
        molecular_weights = gerg2008f.parameters.m
        r = gerg2008f.parameters.r

        x = self.set_concentration(composition)
        print(x)

        m = thermo_props(x, molecular_weights)

        density_r, temperature_r = gerg2008f.reducing_funcs(x)

        delta = density / density_r
        tau = temperature_r / temperature

        ao = gerg2008f.ideal_term(
            x, density, temperature, density_r, temperature_r
        )
        ar = gerg2008f.residual_term(x, delta, tau)

        z = thermo_props.zeta(delta, ar)
        cv = thermo_props.isochoric_heat(tau, r, ao, ar)
        cp = thermo_props.isobaric_heat(delta, tau, r, ao, ar)
        w = thermo_props.sound_speed(delta, tau, r, temperature, m, ao, ar)
        isothermal_thermal_coefficent = (
            thermo_props.isothermal_thermal_coefficent(delta, tau, density, ar)
        )
        dp_dt = thermo_props.dp_dt(density, delta, tau, r, ar)
        dp_drho = thermo_props.dp_drho(temperature, delta, r, ar)
        dp_dv = thermo_props.dp_dv(density, delta, temperature, r, ar)
        p = thermo_props.pressure(delta, density, r, temperature, ar)
        s = thermo_props.entropy(tau, r, ao, ar)
        u = thermo_props.internal_energy(tau, r, temperature, ao, ar)
        h = thermo_props.enthalpy(delta, tau, r, temperature, ao, ar)
        g = thermo_props.gibbs_free_energy(delta, r, temperature, ao, ar)
        jt = thermo_props.joule_thomson_coeff(delta, tau, density, r, ao, ar)
        k = thermo_props.isentropic_exponent(delta, tau, ao, ar)
        b = thermo_props.second_thermal_virial_coeff(density_r, ar)
        c = thermo_props.third_thermal_virial_coeff(density_r, ar)

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
