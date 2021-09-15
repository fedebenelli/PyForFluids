"""PyForFluids
   ===========

   Module that provides.....
"""
import matplotlib.pyplot as plt
import gerg2008


class Model:
    """Thermodynamical model to be applied"""

    def __init__(self, model):
        """Initialize the model
        arguments:
            model (string): Model to be used
                - Possible models: GERG
        """
        self.model = model

        possible_models = {"GERG"}

        if model not in possible_models:
            raise ValueError

        if model == "GERG":
            self.vcalc = gerg2008.vcalc
            self.pressure = gerg2008.thermo_props.pressure
            self.enthalpy = gerg2008.thermo_props.enthalpy
            self.reducing_funcs = gerg2008.reducing_funcs

    def validate_components(self, components):
        """Validate if the components are available
        in the model"""

        non_restricted_models = {"PR", "RKPR"}

        if self.model in non_restricted_models:
            return

        possible_components = {"GERG": ["methane", "carbon_dioxide", ...]}

        for component in components:
            if component not in possible_components[self.model]:
                raise ValueError

    def __repr__(self):
        return "__repr__"


class Fluid(Model):
    """Docstring for GERG."""

    def __init__(self, temperature, pressure=None, density=None, **kwargs):
        """To be defined."""
        if (pressure, density) == (None, None) or None not in (pressure, density):
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
            density = self.vcalc(concentration, pressure, temperature)
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
        self._enthalpy = self.enthalpy(
            self._concentration, self._density, self._temperature
        )

    def isotherm_pv(self, P_range):
        V = self.vcalc(self._concentration, P_range, self._temperature)
        plt.plot(V, P_range)


fluid = Fluid(model="GERG", pressure=250, temperature=150, methane=1, normalize=True)
