"""
"""
from abc import abstractmethod
from functools import partial

import jax.numpy as np
from jax import jit

from pyforfluids.constants import R
from pyforfluids.models.core_models import ArModel, CubicMixingRule


class CubicEoS(ArModel):
    """Generic CubicEoS.

    Base class that contains the logic to calculate the residual Helmhotlz
    energy with a generic Cubic Equation of State.
    """

    autodiff = True

    def set_mixingrule(self, mixingrule: CubicMixingRule):
        self.mixing_rule = mixingrule

    def initial_volume(self, z, pressure, temperature, root="gas"):
        """Initial volume specifier.

        Determine the initial volume to use on a volume solver, based the
        desired root to obtain.

        Parameters
        ----------
        z: array
            Concentrations vector
        pressure: float
            Pressure where to calculate the volume
        temperature:
            Temperature where to calculate the volume
        root: string
            Desired root to obtain, defaults to gas
        """
        if root == "gas":
            return R * temperature / pressure
        elif root == "liquid":
            b = self.repulsive_parameter(z, 1, temperature)
            v0 = self.mixing_rule.mix_b(z, 1, temperature, b)
            return v0

    @partial(jit, static_argnames=["self"])
    def residual_helmholtz(self, z, volume, temperature):
        """Residual Helmholtz method.

        Calculate the residual Helmholtz energy based on the generic cubic
        equation of state:
        """
        v = volume
        t = temperature

        # Pure compounds parameters
        apures = self.attractive_parameter(z, v, t)
        bpures = self.repulsive_parameter(z, v, t)
        delta1 = self.delta1_parameter(z, v, t)
        delta2 = self.delta2_parameter(z, v, t)

        # Mixture "as a fluid" parameters
        a = self.mixing_rule.mix_a(z, v, t, apures)
        b = self.mixing_rule.mix_b(z, v, t, bpures)
        del1 = self.mixing_rule.mix_delta1(z, v, t, delta1)
        del2 = self.mixing_rule.mix_delta2(z, v, t, delta2)

        n = z.sum()
        b_v = b / v

        residual_helmholtz = (
            -np.log(1 - b_v)
            - a/(n*R*t*b) * 1 / (del1 - del2)
            * np.log((1 + del1*b_v)/(1 + del2*b_v))
        )
        # R*t

        return residual_helmholtz

    @abstractmethod
    def attractive_parameter(self, z, v, t):
        raise NotImplementedError

    @abstractmethod
    def repulsive_parameter(self, z, v, t):
        raise NotImplementedError

    @abstractmethod
    def delta1_parameter(self, z, v, t):
        raise NotImplementedError

    @abstractmethod
    def delta2_parameter(self, z, v, t):
        raise NotImplementedError

    @abstractmethod
    def volume_shift(self, z, v, t):
        raise NotImplementedError
