"""
"""
from jax import grad, jacrev, jit

from pyforfluids.constants import R

from scipy.optimize import root_scalar


def get_grad(fun, nvar):
    return jit(grad(fun, argnums=nvar))


def get_hess(fun, nvar_1, nvar_2):
    return jit(jacrev(jit(grad(fun, argnums=nvar_1)), argnums=nvar_2))


class System:
    def __init__(self, armodel=None, gemodel=None, idealmodel=None):
        self.armodel = armodel
        self.gemodel = gemodel
        self.idealmodel = idealmodel

        if armodel:
            self.__setup_armodel_functions(armodel)

    def __setup_armodel_functions(self, armodel):
        if armodel.autodiff:
            # Extract the Ar function
            ar_fun = armodel.residual_helmholtz

            # Function valuation
            self.ar = jit(ar_fun)

            # Only compositional derivatives
            self.dar_dn = get_grad(ar_fun, 0)

            # Only volume derivative
            self.dar_dv = get_grad(ar_fun, 1)

            # Only temperature derivative
            self.dar_dt = get_grad(ar_fun, 2)

            # Second derivatives
            self.dar2_dn2 = get_hess(ar_fun, 0, 0)
            self.dar2_dtn = get_hess(ar_fun, 1, 0)
            self.dar2_dvn = get_hess(ar_fun, 2, 0)

            self.dar2_dt2 = get_hess(ar_fun, 1, 1)
            self.dar2_dtv = get_hess(ar_fun, 2, 1)
            self.dar2_dv2 = get_hess(ar_fun, 2, 2)

    def ln_phi(self, z, volume, temperature):
        ...

    def pressure(self, z, volume, temperature):
        dar_dv = self.dar_dv(z, volume, temperature)
        pressure = -dar_dv + z.sum() * R * temperature / volume
        return pressure

    def volume(self, z, pressure, temperature, root="gas", v0=None):
        def obj_pressure(volume):
            dar_dv = self.dar_dv(z, volume, temperature)
            dar2_dv2 = self.dar2_dv2(z, volume, temperature)

            rt = R * temperature
            rt_v = rt / volume
            n = z.sum()

            p = -dar_dv + n * rt_v
            dp = -rt * dar2_dv2 - n * rt_v / volume

            return p - pressure, dp / 100

        if not v0:
            v0 = self.armodel.initial_volume(z, pressure, temperature, root)

        volume = root_scalar(
            obj_pressure,
            method="newton",
            fprime=True,
            x0=v0,
            xtol=1e-6,
            maxiter=100,
        ).root

        return volume
