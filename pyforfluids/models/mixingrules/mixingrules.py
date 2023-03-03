from functools import partial
from jax import jit
import jax.numpy as np


class ClassicVdW:

    def __init__(self, kij, lij):
        self._kij = kij
        self._lij = lij

    @partial(jit, static_argnames=["self"])
    def mix_a(self, z, v, t, a):
        n = len(a)
        amix = 0

        for i in range(n):
            for j in range(n):
                amix += (
                        z[i] * z[j] 
                        * np.sqrt(a[i] * a[j]) * (1 - self._kij[i, j])
                    )
            
        return amix

    @partial(jit, static_argnames=["self"])
    def mix_b(self, z, v, t, b):
        n = len(b)
        bmix = 0

        for i in range(n):
            for j in range(n):
                bmix += z[i] * z[j] * (b[i] + b[j])/2 * (1 - self._lij[i, j])

        bmix = bmix/np.sum(z)
        
        return bmix

    @partial(jit, static_argnames=["self"])
    def mix_c(self, z, v, t, c):
        return np.sum(z*c)
