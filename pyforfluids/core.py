"""PyForFluids
   ===========

   Module that does stuff.....
"""
import numpy as np


class Fluid:
    def __init__(
        self, model, composition, temperature, pressure=None, density=None
    ):
        if (pressure, density) == (None, None) or None not in (
            pressure,
            density,
        ):
            raise ValueError("Presión y densidad no pueden ser ambos None....")

        self.composition = composition
        self.model = model
        self.temperature = temperature
        self.pressure = pressure
        self.density = density
        self.properties = {}

        model.validate_components(composition)

        self.calculate_properties()

    def set_composition(self, composition):
        """Change the fluid's composition

        Parameter
        ---------
        composition: dict
            Dictionary with the fluid compounds as keys and their molar
            concentration as values

            example:
                composition = {'methane': 0.1, 'ethane': 0.9}
        """
        self.composition = composition

    def set_temperature(self, temperature):
        """Change the fluid's density

        Parameter
        ---------
        temperature: float
            New temperature
        """
        self.temperature = temperature

    def set_density(self, density):
        """Change the fluid's density

        Parameter
        ---------
        density: float
            New density
        """
        self.density = density

    def set_pressure(self, pressure):
        """Change the fluid's pressure

        Parameter
        ---------
        pressure: float
            New pressure

        """
        self.pressure = pressure

    def calculate_properties(self):
        """Calculate the fluid's properties"""
        self.properties = self.model.calculate_properties(
            self.temperature, self.pressure, self.density, self.composition
        )

    def isotherm(self, density_range):
        """Calculate the properties at the fluid temperature along a
        density range.

        Parameters
        ----------
        density_range: array-like
            Range of values of density where to calculate the properties
        """

        # Define a new temporary fluid with the properties of the original one
        fluid = Fluid(
            self.model,
            self.composition,
            self.temperature,
            self.pressure,
            self.density,
        )

        # Initialize a dictionary that will keep all the properties along
        #  the density range
        isotherm = dict()
        fluid.calculate_properties()
        isotherm["density"] = density_range
        for prop in fluid.properties:
            isotherm[prop] = []

        for density in density_range:
            fluid.set_density(density)
            fluid.calculate_properties()

            for prop in fluid.properties:
                isotherm[prop].append(fluid.properties[prop])

        return isotherm

    def __getitem__(self, key):
        return self.properties[key]

    def __repr__(self):
        rep = ""
        rep += "Fluid\n\n"
        rep += "----------------\n"
        rep += "Composition\n"

        for compound in self.composition:
            rep += f"{compound}\t:{self.composition[compound]}\n"
        rep += "----------------\n\n"

        for prop in self.properties:
            prop_value = np.round(self.properties[prop], 6)
            rep += f"{prop}\t: {prop_value}\n"

        return rep
