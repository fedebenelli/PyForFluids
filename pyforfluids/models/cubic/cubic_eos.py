"""
"""
from functools import partial

import jax.numpy as np
from jax import jit

from pyforfluids.constants import R
from pyforfluids.models.residual_helmholtz import ArModel


class CubicEoS(metaclass=ArModel):
    """Generic CubicEoS"""

    autodiff = True

    def __init__(self, equation, mixing_rule):
        self.equation = equation
        self.mixing_rule = mixing_rule

    @partial(jit, static_argnames=["self"])
    def residual_helmholtz(self, z, v, t):
        apures = self.equation.attractive_parameter(z, v, t)
        bpures = self.equation.repulsive_parameter(z, v, t)
        # cpures = self.equation.volume_shift(z, v, t)

        a = self.mixing_rule.mix_a(z, v, t, apures)
        b = self.mixing_rule.mix_b(z, v, t, bpures)
        # c = self.mixing_rule.mix_c(z, v, t, cpures)

        del1 = self.armodel._del1[0]
        del2 = self.armodel._del2[0]

        residual_helmholtz = (
            (
                -np.sum(z) * np.log(1.0 - b / v)
                - a
                / (R * t * b)
                * 1.0
                / (del1 - del2)
                * np.log((1 + del1 * b / v) / (1 + del2 * b / v))
            )
            * R
            * t
        )

        return residual_helmholtz
