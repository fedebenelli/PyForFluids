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
		raise NotImplementhedError()

	@abc.abstractmethod
	def v_calc(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def a_oio(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def a_oir(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def a_ijr(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def ideal_term(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def residual_term(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def pressure(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def enthalpy(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def zeta(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def isochoric_heat(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def isobaric_heat(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def sound_speed(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def isothermal_thermal_coefficent(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def dp_dt(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def dp_drho(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def dp_dv(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def pressure(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def entropy(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def internal_energy(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def enthalpy(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def gibbs_free_energy(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def joule_thomson_coeff(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def isentropic_exponent(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def second_thermal_virial_coeff(self):
		raise NotImplementhedError()

	@abc.abstractmethod
	def third_thermal_virial_coeff(self):
		raise NotImplementhedError()	
