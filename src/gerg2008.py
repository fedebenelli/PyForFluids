from . import core
from .fortran import gerg as gerg2008f
from .fortran import thermoprops as thermopropsf


class GERG2008(core.Model):

    name = "GERG2008"

    # gerg file
    def reducing_funcs(self, x):
        density_r, temperature_r = gerg2008f.reducing_funcs(x)
        return density_r, temperature_r

    def v_calc(self, x, p, temperature):
        density = gerg2008f.V_calc(x, p, temperature)
        return density

    def a_oio(self, density, temperature, density_c, temperature_c, n, v):
        aoio = gerg2008f.a_oio(density, temperature, density_c, temperature_c, n, v)
        return aoio

    def a_oir(self, delta, tau, kpol, kexp, n, d, t, c):
        a_oir = gerg2008f.a_oir(delta, tau, kpol, kexp, n, d, t, c)
        return a_oir

    def a_ijr(self, delta, tau, kpolij, kexpij, n, d, t, eta, eps, gamm, beta):
        a_ijr = gerg2008f.a_ijr(
            delta, tau, kpolij, kexpij, n, d, t, eta, eps, gamm, beta
        )
        return a_ijr

    def ideal_term(self, x, density, temperature, density_r, temperature_r):
        ao = gerg2008f.ideal_term(x, density, temperature, density_r, temperature_r)
        return ao

    def residual_term(self, x, delta, tau):
        ar = gerg2008f.residual_term(x, delta, tau)
        return ar

    # Thermoprops file
    def zeta(self, ar):
        z = thermopropsf.zeta(ar)
        return z

    def isochoric_heat(self, tau, r, ao, ar):
        cv = thermopropsf.isochoric_heat(tau, r, ao, ar)
        return cv

    def isobaric_heat(self, delta, tau, r, ao, ar):
        cp = thermopropsf.isobaric_heat(delta, tau, r, ao, ar)
        return cp

    def sound_speed(self, delta, tau, r, temperature, m, ao, ar):
        w = thermopropsf.sound_speed(delta, tau, r, temperature, m, ao, ar)
        return w

    def isothermal_thermal_coefficent(self, delta, tau, density, ar):
        delta_t = thermopropsf.isothermal_thermal_coefficent(
            delta, tau, density, ar
        )
        return delta_t

    def dp_dt(self, density, delta, tau, r, ar):
        dp_dt = thermopropsf.dp_dt(density, delta, tau, r, ar)
        return dp_dt

    def dp_drho(self, temperature, delta, r, ar):
        dp_drho = thermopropsf.dp_drho(temperature, delta, r, ar)
        return dp_drho

    def dp_dv(self, density, delta, temperature, r, ar):
        dp_dv = thermopropsf.dp_dv(density, delta, temperature, r, ar)
        return dp_dv

    def pressure(self, delta, density, r, temperature, ar):
        p = thermopropsf.pressure(delta, density, r, temperature, ar)
        return p

    def entropy(self, tau, r, ao, ar):
        s = thermopropsf.entropy(tau, r, ao, ar)
        return s

    def internal_energy(self, tau, r, temperature, ao, ar):
        u = thermopropsf.internal_energy(tau, r, temperature, ao, ar)
        return u

    def enthalpy(self, delta, tau, r, temperature, ao, ar):
        h = thermopropsf.enthalpy(delta, tau, r, temperature, ao, ar)
        return h

    def gibbs_free_energy(self, delta, r, temperature, ao, ar):
        g = thermopropsf.gibbs_free_energy(delta, r, temperature, ao, ar)
        return g

    def joule_thomson_coeff(self, delta, tau, density, r, ao, ar):
        jt = thermopropsf.joule_thomson_coeff(delta, tau, density, r, ao, ar)
        return jt

    def isentropic_exponent(self, delta, tau, ao, ar):
        k = thermopropsf.isentropic_exponent(delta, tau, ao, ar)
        return k

    def second_thermal_virial_coeff(self, density_r, ar):
        b = thermopropsf.second_thermal_virial_coeff(density_r, ar)
        return b

    def third_thermal_virial_coeff(self, density_r, ar):
        c = thermopropsf.third_thermal_virial_coeff(density_r, ar)
        return c
