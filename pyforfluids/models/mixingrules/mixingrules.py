from functools import partial

import jax.numpy as np
from jax import jit


class ClassicVdW:
    """ """

    def __init__(self, kij, lij):
        self._kij = np.array(kij)
        self._lij = np.array(lij)

    @partial(jit, static_argnames=["self"])
    def mix(self, z, v, t, a, b, c, del1, del2):
        n = len(z)
        amix = 0
        bmix = 0
        cmix = 0
        del1mix = 0
        del2mix = 0

        for i in range(n):
            for j in range(i, n):
                amix += (
                    2
                    * z[i]
                    * z[j]
                    * np.sqrt(a[i] * a[j])
                    * (1 - self._kij[i, j])
                )

                # Not dividing by two since it's multiplied by two implicitly
                bmix += z[i] * z[j] * (b[i] + b[j]) * (1 - self._lij[i, j])

            cmix += z[i] * c[i]

        del1mix = del1[0]
        del2mix = del1[0]
        bmix = bmix / np.sum(z)

        return amix, bmix, cmix, del1mix, del2mix

    @partial(jit, static_argnames=["self"])
    def mix_a(self, z, v, t, a):
        n = len(z)
        amix = 0

        for i in range(n):
            for j in range(n):
                amix += (
                    z[i] * z[j] * np.sqrt(a[i] * a[j]) * (1 - self._kij[i, j])
                )

        return amix

    @partial(jit, static_argnames=["self"])
    def mix_b(self, z, v, t, b):
        n = len(b)
        bmix = 0

        for i in range(n):
            for j in range(n):
                bmix += z[i] * z[j] * (b[i] + b[j]) / 2 * (1 - self._lij[i, j])

        bmix = bmix / np.sum(z)

        return bmix

    @partial(jit, static_argnames=["self"])
    def mix_c(self, z, v, t, c):
        return np.sum(z * c)

    def mix_delta1(self, z, v, t, delta1):
        return delta1[0]

    def mix_delta2(self, z, v, t, delta2):
        return delta2[0]
