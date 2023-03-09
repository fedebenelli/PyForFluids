from functools import partial

import jax.numpy as np
from jax import jit

from pyforfluids.constants import R

from pyforfluids.models.cubic.cubic_eos import CubicEoS


class PR76(CubicEoS):
    """
    Peng Robinson 76 Equation of State
    """

    def __init__(self, tc, pc, w, mixing_rule):
        self.set_mixingrule(mixing_rule)
        self._del1 = np.array([1 + np.sqrt(2.0)])
        self._del2 = np.array([1 - np.sqrt(2.0)])

        self._tc = np.array(tc)
        self._pc = np.array(pc)
        self._w = np.array(w)

        self.n = len(tc)

        self._ac = 0.45723553 * R**2 * tc**2 / pc
        self._b = 0.07779607 * R * tc / pc
        self._k = 0.37464 + 1.54226 * w - 0.26993 * w**2
        self._c = np.zeros(self.n)

    @partial(jit, static_argnames=["self"])
    def attractive_parameter(self, z, v, t):
        a = self._ac * (1 + self._k * (1 - np.sqrt(t / self._tc))) ** 2
        return a

    @partial(jit, static_argnames=["self"])
    def repulsive_parameter(self, z, v, t):
        return self._b

    @partial(jit, static_argnames=["self"])
    def volume_shift(self, z, v, t):
        return self._c

    @partial(jit, static_argnames=["self"])
    def delta1_parameter(self, z, v, t):
        return self._del1

    @partial(jit, static_argnames=["self"])
    def delta2_parameter(self, z, v, t):
        return self._del2
