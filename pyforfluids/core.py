"""Core module."""
import warnings

import numpy as np


class Fluid:
    """Describes a fluid based on a given model and it's thermo variables.

    Density and pressure can't be defined at the same time. If pressure is
    given, the density will be calculated with an iterative algorithm using
    the derivatives given by the model.

    Parameters
    ----------
    model: pyforfluids model_like
        Model to use in the properties calculation.
    composition : dict
        Dictionary with the compounds concentrations as:
        ``{'methane': 0.8, 'ethane': 0.1}``
        In some cases, as in GERG2008, the values will be normalized for the
        calculations but won't be modified in the Fluid attribute
    temperature: float
        Fluid temperature in degrees Kelvin [K]
    pressure: float
        Fluid pressure in Pascals [Pa]
    density: float
        Fluid density in mol per liter [mol/L]

    Attributes
    ----------
    properties : dict
        Fluid properties calculated by it's model.

    Methods
    -------
    copy:
        Returns a copy of the Fluid.
    set_composition:
        Change the Fluid's composition.
    set_temperature:
        Change the Fluid's temperature.
    set_density:
        Change the Fluid's density.
    set_pressure:
        Change the FLuid's pressure.
    calculate_properties:
        Calculate the Fluids properties, returns as a dictionary.
    isotherm:
        Calculate the Fluid properties along a density range.
    density_iterator:
        Calculate the Fluid's density based on a specified pressure.
    """

    def __init__(
        self, model, composition, temperature, pressure=None, density=None
    ):
        if (pressure, density) == (None, None) or None not in (
            pressure,
            density,
        ):
            raise ValueError(
                "Can't define pressure and density at the same time"
            )

        self.composition = composition
        self.model = model
        self.temperature = temperature
        self.pressure = pressure
        self.density = density
        self.properties = {}

        model.validate_components(composition)

        if density is None:
            # If no density was given calculate it
            self.density = np.nan
            density = self.density_iterator(pressure)[0]
            self.density = density

        self.calculate_properties()

    def copy(self):
        """Return a copy of the fluid, taking density as independant variable.

        Returns
        -------
        Fluid
        """
        return Fluid(
            model=self.model,
            composition=self.composition,
            temperature=self.temperature,
            density=self.density,
        )

    def set_composition(self, composition):
        """Change the fluid's composition.

        Parameters
        ----------
        composition: dict
            Dictionary with the fluid compounds as keys and their molar
            concentration as values

            example:
            composition = {'methane': 0.1, 'ethane': 0.9}
        """
        self.composition = composition

    def set_temperature(self, temperature):
        """Change the fluid's temperature.

        Parameters
        ----------
        temperature: float
            New temperature
        """
        self.temperature = temperature

    def set_density(self, density):
        """Change the fluid's density.

        Parameters
        ----------
        density: float
            New density
        """
        self.density = density

    def set_pressure(self, pressure):
        """Change the fluid's pressure.

        Parameters
        ----------
        pressure: float
            New pressure
        """
        self.pressure = pressure
        new_density = self.density_iterator(pressure)[0]
        self.density = new_density

    def calculate_properties(self):
        """Calculate the fluid's properties."""
        self.properties = self.model.calculate_properties(
            self.temperature, self.pressure, self.density, self.composition
        )

    def isotherm(self, density_range):
        """Calculate isotherm along a density range.

        Calculate the fluid's properties that it's model can give at constant
        temperature along a density range.

        Parameters
        ----------
        density_range: array-like
            Range of values of density where to calculate the properties.

        Returns
        -------
        dict
            Dictionary with all the properties that the model can calculate
            along the density_range.
        """
        # Define a new temporary fluid with the properties of the original one
        fluid = self.copy()

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

    def density_iterator(self, objective_pressure):
        """With a given pressure and temperature value, get the fluid density.

        Parameters
        ----------
        objective_pressure: float
            Fluid pressure where to calculate density.
        temperature: float
            Fluid temperature where to calculate density.

        Returns
        -------
        float
            Calculated density.
        float
            Pressure where the density converged.
        int
            Number of iterations.
        """
        fluid = self.copy()

        t = fluid.temperature

        fluid.set_temperature(t)

        step = 1
        it = 0
        r = 8.31446261815324

        rho_i = objective_pressure / (r * t) / 1000

        fluid.set_density(rho_i)
        fluid.calculate_properties()

        p = fluid.properties["p"]
        precision = 0.01

        while abs(p - objective_pressure) > objective_pressure * precision:
            it = it + 1

            # Calculate properties on the new Newton point
            fluid.set_density(rho_i)
            fluid.calculate_properties()
            p = fluid.properties["p"]
            dp_drho = fluid.properties["dp_drho"] * 1000

            delta = (p - objective_pressure) / dp_drho

            rho_i = rho_i - step * delta

            if it > 50:
                warnings.warn(
                    RuntimeWarning("Couldn't converge with 50 iterations")
                )
                return rho_i, p, it

        return rho_i, p, it

    def __getitem__(self, key):
        """Access the fluid properties as a dictionary."""
        return self.properties[key]

    def __repr__(self):
        """Give a summary table of the fluid properties."""
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
