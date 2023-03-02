from functools import partial

from jax import jit

import jax.numpy as np

from scipy.constants import R


class PR76:
    def __init__(self, tc, pc, w):
        self.n = len(tc)
        self._del1 = [1 + np.sqrt(2.0)]
        self._del2 = [1 - np.sqrt(2.0)]

        self._tc = tc
        self._pc = pc
        self._w = w

        self._ac = 0.45723553 * R**2 * tc**2 / pc
        self._b = 0.07779607*R*tc/pc
        self._k = 0.37464 + 1.54226*w - 0.26993*w**2
        self._c = 0

    @partial(jit, static_argnames=["self"])
    def attractive_parameter(self, z, v, t):
        a = self._ac * (1 + self._k*(1 - np.sqrt(t/self._tc)))**2
        return a
    
    @partial(jit, static_argnames=["self"])
    def repulsive_parameter(self, z, v, t):
        return self._b
    
    @partial(jit, static_argnames=["self"])
    def volume_shift(self, z, v, t):
        return self._c
    
    @partial(jit, static_argnames=["self"])
    def del1(self, z, v, t):
        return self._del1

    @partial(jit, static_argnames=["self"])
    def del2(self, z, v, t):
        return self._del2
