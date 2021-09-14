"""DOCSTRING
"""
import gerg2008


class GERG(object):
    """Docstring for GERG. """

    def __init__(self, temperature, pressure=None, density=None, **kwargs):
        """TODO: to be defined.


        :temperature: TODO
        :pressure: TODO
        :density: TODO
        :**kwargs: TODO

        """
        if (pressure, density) == (None,None) or None not in (pressure, density):
            raise ValueError("Presión y densidad no pueden ser ambos None....")

        self.methane = kwargs.get("methane", 0)
        self.nitrogen = 0
        self.carbon_dioxide = 0
        self.ethane = 0
        self.propane = 0
        self.butane = 0
        self.isobutane = 0
        self.pentane = 0
        self.isopentane = 0
        self.hexane = 0
        self.heptane = 0
        self.octane = 0
        self.nonane = 0
        self.decane = 0
        self.hydrogen = 0
        self.oxygen = 0
        self.carbon_monoxide = 0
        self.water = 0
        self.hydrogen_sulfide = 0
        self.helium = 0
        self.argon = 0

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

        # Convertir unidades a SI
        ...

        if density is None:
            # Get density with concentration, pressure and temperature
            density = gerg2008.vcalc(concentration, pressure, temperature)
        else:
            pressure = gerg2008.pressure(concentration, density, temperature)

        self._concentration = concentration
        self._temperature = temperature
        self._pressure = pressure
        self._density = density
        self._enthalpy = 0
        self._entropy = 0

    def get_properties(self):
        self._enthalpy = gerg2008.enthalpy(self.concentration, self.density, self.temperature)
        self._entropy = gerg2008.entropy(self.concentration, density, temperature)





