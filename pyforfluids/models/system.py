"""
"""
from functools import partial

from jax import jit

# import numpy as np
import jax.numpy as np

from pyforfluids.autodiff import get_grad, get_hess
from pyforfluids.constants import R
from pyforfluids.models.core_models import ArModel, GeModel

from scipy.optimize import newton


class System:
    """Thermodynamic System

    Object that orchestrates different kind of thermodynamic models defined
    in PyForFluids.

    A system is assumed to contain one or more of the following attributes:

    - `ArModel`: Residual helmholtz model that calculates the residual Helmholtz
       energy and all it's relevant derivatives.
    - `GeModel`: Residual helmholtz model that calculates the Excess Gibbs
        energy and all it's relevant derivatives.
    - `IdealModel`: Residual helmholtz model that calculates the ideal Helmholtz
        energy and all it's relevant derivatives.

    Parameters
    ----------
    armodel: ArModel
        Object that calculates the residual Helmholtz energy and it's
        derivatives.
    gemodel: GeModel
        Object that calculates the excess Gibbs energy and it's
        derivatives.
    idealmodel: IdealModel
        Object that calculates the ideal Helmholtz energy and it's
        derivatives.

    All the models can make the use of automatic differentiation but if the
    analytical derivatives are provided they can be used as long as the
    general interface is followed.

    The thermodynamic identities are calculated as defined by the book
    of Michelsen and Mollerup [cite].
    """

    def __init__(
        self, armodel: ArModel = None, gemodel: GeModel = None, idealmodel=None
    ):
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

    def compressibility_factor(self, z, volume, temperature):
        pressure = self.pressure(z, volume, temperature)
        return pressure * volume / (R * temperature)

    def ln_phi(self, z, volume, temperature):
        """Fugacity coefficent"""
        z_fact = self.compressibility_factor(z, volume, temperature)
        dar_dn = self.dar_dn(z, volume, temperature)
        lnphi = dar_dn - np.log(z_fact)

        return lnphi

    def pressure(self, z, volume, temperature):
        """Pressure

        Obtain the pressure as a function of composition, volume and temperature

        Parameters
        ----------
        z: array
            Vector of compositions
        volume: float
            Volume [L]
        temperature: float
            temperature [K]
        """
        moles = z.sum()
        v = volume
        t = temperature

        dar_dv = self.dar_dv(z, v, t)

        p = R * t * (moles / v - dar_dv)

        return p

    def dp_dv(self, z, volume, temperature):
        moles = z.sum()
        v = volume
        t = temperature

        dar2_dv2 = self.dar_dv(z, v, t)

        rt = R * t
        dpdv = -rt * dar2_dv2 - moles * rt / v**2

        return dpdv

    def volume(
        self, z, pressure, temperature, root="gas", v0=None, tol=1e-3, it=50
    ):
        """Volume solver as a function of pressure and temperature"""

        if not v0:
            # Find the initial volume using the model's initializer
            v0 = self.armodel.initial_volume(
                z, pressure, temperature, root=root
            )

        v = v0 * 1.1

        vmin = v0
        vmax = v0 * 100
        delta = 1

        i = 1

        for i in range(1, it + 1):
            p = self.pressure(z, v, temperature)

            # Calculate Newton step
            dp = self.dp_dv(z, v, temperature)
            delta = (pressure - p) / dp

            if root == "liquid":
                # Corrections to find the liquid root
                if p > pressure:
                    vmin = v
                else:
                    vmax = v

                while abs(delta) > 0.9 * v:
                    delta /= 2

            v = v + delta

            if v > vmax or v < vmin and root == "liquid":
                v = (vmin + vmax)/2

            if abs(p - pressure) < tol:
                break
        return v, i
