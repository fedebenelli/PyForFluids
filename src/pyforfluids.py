"""PyForFluids
   ===========

   Module that does stuff.....
"""

import matplotlib.pyplot as plt

import gerg2008


class Model:
    """Determine the model to be utilized"""

    def __init__(self, model):
        """Initialize the model
        arguments:
            model (string): Model to be used
                - Possible models: GERG
        """
        self.model = model

        possible_models = {"GERG"}

        if model not in possible_models:
            raise ValueError(f"{model} isn't implemented, yet")

        if model == "GERG":
            # Define the model's methods
            self.reducing_funcs = gerg2008.reducing_funcs
            self.v_calc = gerg2008.v_calc
            self.a_oio = gerg2008.a_oio
            self.a_oir = gerg2008.a_oir
            self.a_ijr = gerg2008.a_ijr
            self.ideal_term = gerg2008.ideal_term
            self.residual_term = gerg2008.residual_term

            self.pressure = gerg2008.thermo_props.pressure
            self.enthalpy = gerg2008.thermo_props.enthalpy
            self.zeta = gerg2008.thermo_props.zeta
            self.isochoric_heat = gerg2008.thermo_props.isochoric_heat
            self.isobaric_heat = gerg2008.thermo_props.isobaric_heat
            self.sound_speed = gerg2008.thermo_props.sound_speed
            self.isothermal_thermal_coefficent = (
                gerg2008.thermo_props.isothermal_thermal_coefficent
            )
            self.dp_dt = gerg2008.thermo_props.dp_dt
            self.dp_drho = gerg2008.thermo_props.dp_drho
            self.dp_dv = gerg2008.thermo_props.dp_dv
            self.pressure = gerg2008.thermo_props.pressure
            self.entropy = gerg2008.thermo_props.entropy
            self.internal_energy = gerg2008.thermo_props.internal_energy
            self.enthalpy = gerg2008.thermo_props.enthalpy
            self.gibbs_free_energy = gerg2008.thermo_props.gibbs_free_energy
            self.joule_thomson_coeff = (
                gerg2008.thermo_props.joule_thomson_coeff
            )
            self.isentropic_exponent = (
                gerg2008.thermo_props.isentropic_exponent
            )
            self.second_thermal_virial_coeff = (
                gerg2008.thermo_props.second_thermal_virial_coeff
            )
            self.third_thermal_virial_coeff = (
                gerg2008.thermo_props.third_thermal_virial_coeff
            )

    def validate_components(self, components):
        """Validate if the components are available
        in the model"""

        non_restricted_models = {"PR", "RKPR"}

        if self.model in non_restricted_models:
            return

        possible_components = {
            "GERG": [
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
        }

        for component in components:
            if component not in possible_components[self.model]:
                raise ValueError(
                    f"{component} Can't be used in model {self.model}"
                )


class Fluid(Model):
    """Docstring for GERG."""

    def __init__(self, temperature, pressure=None, density=None, **kwargs):
        """To be defined."""
        if (pressure, density) == (None, None) or None not in (
            pressure,
            density,
        ):
            raise ValueError("Presión y densidad no pueden ser ambos None....")

        self.methane = kwargs.get("methane", 0)
        self.nitrogen = kwargs.get("nitrogen", 0)
        self.carbon_dioxide = kwargs.get("carbon_dioxide", 0)
        self.ethane = kwargs.get("ethane", 0)
        self.propane = kwargs.get("propane", 0)
        self.butane = kwargs.get("butane", 0)
        self.isobutane = kwargs.get("isobutane", 0)
        self.pentane = kwargs.get("pentane", 0)
        self.isopentane = kwargs.get("isopentane", 0)
        self.hexane = kwargs.get("hexane", 0)
        self.heptane = kwargs.get("heptane", 0)
        self.octane = kwargs.get("octane", 0)
        self.nonane = kwargs.get("nonane", 0)
        self.decane = kwargs.get("decane", 0)
        self.hydrogen = kwargs.get("hydrogen", 0)
        self.oxygen = kwargs.get("oxygen", 0)
        self.carbon_monoxide = kwargs.get("carbon_monoxide", 0)
        self.water = kwargs.get("water", 0)
        self.hydrogen_sulfide = kwargs.get("hydrogen_sulfide", 0)
        self.helium = kwargs.get("helium", 0)
        self.argon = kwargs.get("argon", 0)

        # Definir concentraciones dadas
        #   considerar definir concentraciones como un clase aparte
        concentration = [
            self.methane,
            self.nitrogen,
            self.carbon_dioxide,
            self.ethane,
            self.propane,
            self.butane,
            self.isobutane,
            self.pentane,
            self.isopentane,
            self.hexane,
            self.heptane,
            self.octane,
            self.nonane,
            self.decane,
            self.hydrogen,
            self.oxygen,
            self.carbon_monoxide,
            self.water,
            self.hydrogen_sulfide,
            self.helium,
            self.argon,
        ]

        # Validar que los componentes pertenezcan al modelo
        self.validate_components(kwargs)

        # Convertir unidades a SI
        ...

        if density is None:
            # Get density with concentration, pressure and temperature
            density = self.v_calc(concentration, pressure, temperature)
        else:
            pressure = self.pressure(concentration, density, temperature)

        self._delta, self._tau = self.reducing_funcs(concentration)

        self._concentration = concentration
        self._temperature = temperature
        self._pressure = pressure
        self._density = density
        self._enthalpy = self.enthalpy(concentration, density, temperature)
        self._entropy = 0

    def update_properties(self):
        """ """
        self._enthalpy = self.enthalpy(
            self._concentration, self._density, self._temperature
        )

    def isotherm_pv(self, p_range):
        """ """
        v_range = self.v_calc(self._concentration, p_range, self._temperature)
        plt.plot(v_range, p_range)


fluid = Fluid(
    model="GERG", pressure=250, temperature=150, methane=1, normalize=True
)
