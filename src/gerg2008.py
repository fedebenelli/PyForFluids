from . import core
from .fortran import gerg as gerg2008F
from .fortran import thermoprops as thermopropsF


class GERG2008(core.Model):
	
	name = "GERG2008"
	
	# gerg file
	def reducing_funcs(self, X):
		rho_r = gerg2008F.reducing_funcs(X)[0]
		T_r = gerg2008F.reducing_funcs(X)[1]
		return rho_r, T_r

	
	def v_calc(self,X, P, T):
		rho = gerg2008F.V_calc(X)
		return rho

	
	def a_oio(self,rho, T, rho_c, T_c, n, v):
		aoio = gerg2008F.a_oio(rho, T, rho_c, T_c, n, v)
		return aoio
		

	def a_oir(self,delta, tau, Kpol, Kexp, n, d, t, c):
		a_oir = gerg2008F.a_oir(delta, tau, Kpol, Kexp, n, d, t, c)
		return a_oir

	
	def a_ijr(self, delta, tau, Kpolij, Kexpij, n, d, t, eta, eps, gamm, beta):
		a_ijr = gerg2008F.a_ijr(delta, tau, Kpolij, Kexpij, n, d, t, eta, eps, gamm, beta)
		return a_ijr 

	
	def ideal_term(self,X, rho, T, rho_r, T_r):
		ao = gerg2008F.ideal_term(X, rho, T, rho_r, T_r)
		return ao

	
	def residual_term(self,X, delta, tau):
		ar = gerg2008F.residual_term(X, delta, tau)
		return ar

	# Thermoprops file
	def pressure(self,delta, rho, R, T, Ar):
		p = thermopropsF.pressure(delta, rho, R, T, Ar)
		return p

	
	def enthalpy(self,delta, tau, R, T, Ao, Ar):
		h = thermopropsF.enthalpy(delta, tau, R, T, Ao, Ar)
		return h

	
	def zeta(self, ar):
		z = thermopropsF.zeta(ar)
		return z

	
	def isochoric_heat(self,tau, R, Ao, Ar):
		cv = thermopropsF.isochoric_heat(tau, R, Ao, Ar)
		return cv

	
	def isobaric_heat(self, delta, tau, R, Ao, Ar):
		cp = thermopropsF.isobaric_heat(delta, tau, R, Ao, Ar)
		return cp

	
	def sound_speed(self, delta, tau, R, T, M, Ao, Ar):
		w = thermopropsF.sound_speed(delta, tau, R, T, M, Ao, Ar)
		return w

	
	def isothermal_thermal_coefficent(self,delta, tau, rho, Ar):
		delta_t = thermopropsF.isothermal_thermal_coefficent(delta, tau, rho, Ar)
		return delta_t

	
	def dp_dt(self,rho, delta, tau, R, Ar):
		dp_dt = thermopropsF.dp_dt(rho, delta, tau, R, Ar)
		return dp_dt

	
	def dp_drho(self,T, delta, R, Ar):
		dp_drho = thermopropsF.dp_drho(T, delta, R, Ar)
		return dp_drho

	
	def dp_dv(self,rho, delta, T, R, Ar):
		dp_dv = thermopropsF.dp_dv(rho, delta, T, R, Ar)
		return dp_dv

	
	def pressure(self,delta, rho, R, T, Ar):
		p = thermopropsF.pressure(delta, rho, R, T, Ar)
		return p

	
	def entropy(self,tau, R, Ao, Ar):
		s = thermopropsF.entropy(tau, R, Ao, Ar)
		return s

	
	def internal_energy(self,tau, R, T, Ao, Ar):
		u = thermopropsF.internal_energy(tau, R, T, Ao, Ar)
		return u

	
	def enthalpy(self,delta, tau, R, T, Ao, Ar):
		h = thermopropsF.enthalpy(delta, tau, R, T, Ao, Ar)
		return h

	
	def gibbs_free_energy(self,delta, R, T, Ao, Ar):
		g = thermopropsF.gibbs_free_energy(delta, R, T, Ao, Ar)
		return g
	
	def joule_thomson_coeff(self,delta, tau, rho, R, Ao, Ar):
		JT = thermopropsF.joule_thomson_coeff(delta, tau, rho, R, Ao, Ar)
		return JT

	
	def isentropic_exponent(self,delta, tau, Ao, Ar):
		k = thermopropsF.isentropic_exponent(delta, tau, Ao, Ar)
		return k

	
	def second_thermal_virial_coeff(self,rho_r, Ar):
		B = thermopropsF.second_thermal_virial_coeff(rho_r, Ar)
		return B

	
	def third_thermal_virial_coeff(self,rho_r, Ar):
		C = thermopropsF.third_thermal_virial_coeff(rho_r, Ar)
		return C
	

