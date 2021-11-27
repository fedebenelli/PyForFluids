"""GERG2008 EoS."""
import warnings

import numpy as np

from ..fortran import gerg2008f
from ..fortran.thermo_props import thermo_props


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
        Check in which range of validity are the pressure and density.
    set_concentration:
        Normalize the composition as molar fractios.
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
            raise ValueError(f"'{diff}' ain't valid components")

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
        elif pressure > 70.0:
            warnings.warn(
                "Working conditions belong to the invalid vality range.",
                category=UserWarning,
            )
        elif (temperature >= 90.0 and temperature <= 450.0) and (
            pressure <= 35.0
        ):
            return
        else:
            warnings.warn(
                "Working conditions belong to the extended vality range.",
                category=UserWarning,
            )
            return

    def set_concentration(self, composition):
        """Verify if the sum of the molar fractions of the fluid components is 1.

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
            warnings.warn(
                "Composition doesn't add '1', will be normalized",
                category=UserWarning,
            )
            x = x / sum_value

        return x

    def calculate_properties(
        self, temperature, pressure, density, composition
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
        gerg2008f.get_params()
        molecular_weights = gerg2008f.parameters.m
        r = gerg2008f.parameters.r

        # Concentration dependant parameters
        x = self.set_concentration(composition)
        m = thermo_props.mean_molecular_weight(x, molecular_weights)

        # Reducing functions
        density_r, temperature_r = gerg2008f.reducing_funcs(x)
        delta = density / density_r
        tau = temperature_r / temperature

        # Model tems
        ao = gerg2008f.ideal_term(
            x, density, temperature, density_r, temperature_r
        )
        ar = gerg2008f.residual_term(x, delta, tau)

        # Properties
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
        ar_virial = gerg2008f.residual_term(x, 1e-15, tau)
        b = thermo_props.second_thermal_virial_coeff(density_r, ar_virial)
        c = thermo_props.third_thermal_virial_coeff(density_r, ar_virial)

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
