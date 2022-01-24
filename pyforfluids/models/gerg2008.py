""""""
import warnings

import numpy as np

from ..fortran import gerg2008f
from ..fortran.thermo_props import thermo_props as props


class GERG2008:
    """ """

    name = "GERG2008"

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

    def validate_components(self, components):
        given_components = set(components)

        diff = given_components.difference(self.valid_components)
        if len(diff) > 0:
            warnings.warn(
                f"{self.name} Valid Components:\n{self.valid_components}"
            )
            raise ValueError(f"'{diff}' ain't valid components")

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
        self, temperature, pressure, density, composition, ideal=False
    ):

        # General use parameteres
        gerg2008f.get_params()
        molecular_weights = gerg2008f.parameters.m
        r = gerg2008f.parameters.r

        # Concentration dependant parameters
        x = self.set_concentration(composition)
        m = props.mean_molecular_weight(x, molecular_weights)

        # Reducing functions
        reducing_variables = gerg2008f.reducing_funcs(x)
        reducing_density = reducing_variables[0]
        reducing_temperature = reducing_variables[1]
        dtr_dx = reducing_variables[2]
        dvr_dx = reducing_variables[3]

        delta = density / reducing_density
        tau = reducing_temperature / temperature

        # Model tems
        ao = gerg2008f.ideal_term(
            x, density, temperature, reducing_density, reducing_temperature
        )

        (
            ar,
            ar_x,
            ar_dx,
            ar_tx,
            ar_xx,
        ) = gerg2008f.residual_term(x, delta, tau)

        if ideal:
            ar = np.zeros_like(ar)
            ar_x = np.zeros_like(ar_x)
            ar_dx = np.zeros_like(ar_dx)
            ar_tx = np.zeros_like(ar_tx)
            ar_xx = np.zeros_like(ar_xx)

        # Properties
        z = props.zeta(delta, ar)
        p = props.pressure(delta, density, r, temperature, ar)
        s = props.entropy(tau, r, ao, ar)
        u = props.internal_energy(tau, r, temperature, ao, ar)
        h = props.enthalpy(delta, tau, r, temperature, ao, ar)
        g = props.gibbs_free_energy(delta, r, temperature, ao, ar)
        jt = props.joule_thomson_coeff(delta, tau, density, r, ao, ar)
        k = props.isentropic_exponent(delta, tau, ao, ar)
        cv = props.isochoric_heat(tau, r, ao, ar)
        cp = props.isobaric_heat(delta, tau, r, ao, ar)
        w = props.sound_speed(delta, tau, r, temperature, m, ao, ar)
        isothermal_thermal_coefficent = props.isothermal_thermal_coefficent(
            delta, tau, density, ar
        )

        # PVT derivatives
        dp_dt = props.dp_dt(density, delta, tau, r, ar)
        dp_drho = props.dp_drho(temperature, delta, r, ar)
        dp_dv = props.dp_dv(density, delta, temperature, r, ar)

        # Per-mol derivatives
        dar_dn, dadr_dn = props.helmholtz_per_mol(
            x,
            delta,
            tau,
            reducing_density,
            reducing_temperature,
            ar,
            ar_x,
            ar_dx,
            dvr_dx,
            dtr_dx,
        )
        dnar_dn = ar[0, 0] + dar_dn

        fugacity_coefficent = dnar_dn - np.log(z)
        fugacity = x * density * r * temperature * np.exp(dnar_dn)

        ar_virial, *_ = gerg2008f.residual_term(x, 1e-15, tau)
        b = props.second_thermal_virial_coeff(reducing_density, ar_virial)
        c = props.third_thermal_virial_coeff(reducing_density, ar_virial)

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
            "dar_dn": dar_dn,
            "dadr_dn": dadr_dn,
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
            "fugacity_coefficent": fugacity_coefficent,
            "fugacity": fugacity,
        }
