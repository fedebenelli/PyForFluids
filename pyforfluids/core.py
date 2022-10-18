"""Core module."""
import warnings

import numpy as np

import pandas as pd

from scipy.constants import R
from scipy.optimize import root_scalar


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

        # Validate the fluid's components before any calculation is made
        model.validate_components(composition)

        # If no density was given calculate it and use the gas value
        if density is None:
            self.set_pressure(pressure)

        self.calculate_properties()
        self.pressure = self["pressure"]

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
        self.density = np.nan
        liquid_density, vapor_density, single_fase = self.density_iterator(
            pressure
        )

        if not single_fase:
            warnings.warn(
                "Two roots were found! Vapor-phase value will be used",
                category=UserWarning,
            )

        self.density = vapor_density

    def calculate_properties(self):
        """Calculate the fluid's properties."""
        self.properties = self.model.calculate_properties(
            self.temperature,
            self.pressure,
            self.density,
            self.composition,
        )
        # Update the pressure with the new pressure value
        self.pressure = self.properties["pressure"]

        self.properties = pd.Series(self.properties)

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
        isotherm = {}

        fluid.calculate_properties()
        isotherm["density"] = density_range

        properties = fluid.properties.index

        for prop in properties:
            isotherm[prop] = []

        for density in density_range:
            fluid.set_density(density)
            fluid.calculate_properties()

            for prop in properties:
                isotherm[prop].append(fluid.properties[prop])

        isotherm = pd.DataFrame(isotherm)

        return isotherm

    def density_iterator(
        self, objective_pressure, vapor_phase=True, liquid_phase=True
    ):
        """With a given pressure and temperature value, get the fluid density.

        Parameters
        ----------
        objective_pressure: float
            Fluid pressure where to calculate density.
        vapor_phase: bool
            Find vapor phase root.
        liquid_phase: bool
            Find liquid phase root.

        Returns
        -------
        liquid_density: float
            Calculated density in liquid phase.
        vapor_density: float
            Calculated density in vapor phase.
        single_phase: bool
            True if only one root was found.
        """

        def fluid_pressure(density, fluid, obj):
            temp = fluid.copy()
            temp.set_density(density)
            temp.calculate_properties()

            return temp["pressure"] - obj, temp["dp_drho"] * 1000

        def find_root(fluid, x0, objective_pressure):
            root = root_scalar(
                fluid_pressure,
                args=(fluid, objective_pressure),
                x0=x0,
                fprime=True,
                method="newton",
                xtol=1e-4,
            )
            sol = root.root

            return sol

        fluid = self.copy()

        # LIQUID ROOT
        liquid_density = None
        if liquid_phase:
            initial_density = 25
            liquid_density = find_root(
                fluid, initial_density, objective_pressure
            )

        # ---------------------------------------------------------

        # GAS ROOT
        # Use ideal gas density
        vapor_density = None
        if vapor_phase:
            initial_density = (
                objective_pressure / (R * fluid.temperature) / 1000
            )
            vapor_density = find_root(
                fluid, initial_density, objective_pressure
            )

        if vapor_phase and liquid_phase:
            single_phase = np.allclose(liquid_density, vapor_density)
        else:
            single_phase = True

        return liquid_density, vapor_density, single_phase

    def __getitem__(self, key):
        """Access the fluid properties as a dictionary."""
        return self.properties[key]

    def __repr__(self):
        """Fluid's repr."""
        rep = (
            f"Fluid(model={self.model}, temperature={self.temperature}, "
            f"pressure={self.pressure:.4f}, density={self.density:.4f}, "
            f"composition={self.composition})"
        )
        return rep
