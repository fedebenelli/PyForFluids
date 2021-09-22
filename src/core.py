"""PyForFluids
   ===========

   Module that does stuff.....
"""

import abc

import matplotlib.pyplot as plt


class Model(metaclass=abc.ABCMeta):

    name = None

    def __init_subclass__(cls):
        super().__init_subclass__()
        if cls.name is None:
            raise TypeError("class attribute 'name' must be defined")

    @abc.abstractmethod
    def reducing_funcs(self):
        raise NotImplementedError

    @abc.abstractmethod
    def v_calc(self):
        raise NotImplementedError

    @abc.abstractmethod
    def a_oio(self):
        raise NotImplementedError

    @abc.abstractmethod
    def a_oir(self):
        raise NotImplementedError

    @abc.abstractmethod
    def a_ijr(self):
        raise NotImplementedError

    @abc.abstractmethod
    def ideal_term(self):
        raise NotImplementedError

    @abc.abstractmethod
    def residual_term(self):
        raise NotImplementedError

    @abc.abstractmethod
    def zeta(self):
        raise NotImplementedError

    @abc.abstractmethod
    def isochoric_heat(self):
        raise NotImplementedError

    @abc.abstractmethod
    def isobaric_heat(self):
        raise NotImplementedError

    @abc.abstractmethod
    def sound_speed(self):
        raise NotImplementedError

    @abc.abstractmethod
    def isothermal_thermal_coefficent(self):
        raise NotImplementedError

    @abc.abstractmethod
    def dp_dt(self):
        raise NotImplementedError

    @abc.abstractmethod
    def dp_drho(self):
        raise NotImplementedError

    @abc.abstractmethod
    def dp_dv(self):
        raise NotImplementedError

    @abc.abstractmethod
    def pressure(self):
        raise NotImplementedError

    @abc.abstractmethod
    def entropy(self):
        raise NotImplementedError

    @abc.abstractmethod
    def internal_energy(self):
        raise NotImplementedError

    @abc.abstractmethod
    def enthalpy(self):
        raise NotImplementedError

    @abc.abstractmethod
    def gibbs_free_energy(self):
        raise NotImplementedError

    @abc.abstractmethod
    def joule_thomson_coeff(self):
        raise NotImplementedError

    @abc.abstractmethod
    def isentropic_exponent(self):
        raise NotImplementedError

    @abc.abstractmethod
    def second_thermal_virial_coeff(self):
        raise NotImplementedError

    @abc.abstractmethod
    def third_thermal_virial_coeff(self):
        raise NotImplementedError


class Fluid(Model):
    """Docstring for GERG."""

    def __init__(self, temperature, pressure=None, density=None, **kwargs):
        """To be defined."""
        if (pressure, density) == (None, None) or None not in (
            pressure,
            density,
        ):
            raise ValueError("Presión y densidad no pueden ser ambos None....")

        # Validar que los componentes pertenezcan al modelo
        self.validate_components(kwargs)

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
        self._entropy = self.entropy(concentration, density, temperature)

    def update_properties(self):
        """ """
        self._enthalpy = self.enthalpy(
            self._concentration, self._density, self._temperature
        )

    def isotherm_pv(self, p_range):
        """ """
        v_range = self.v_calc(self._concentration, p_range, self._temperature)
        plt.plot(v_range, p_range)
