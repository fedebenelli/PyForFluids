"""
"""
import jax.numpy as np
from jax import jacrev, jit


class System:
    def __init__(self, armodel, gemodel, idealmodel):
        self.armodel = armodel
        self.gemodel = gemodel
        self.idealmodel = idealmodel

        if armodel.autodiff:
            self._residual_helmholtz = jit(jacrev(armodel.residual_helmholtz))

    def residual_helmholtz(self, z, volume, temperature):
        return self._residual_helmholtz(z, volume, temperature)

    def ln_phi(self, z, volume, temperature):
        # residual_helmholtz = self.residual_helmholtz(z, volume, temperature)
        ...

    def pressure(self, z, volume, temperature):
        ...

    def volume(self, z, pressure, temperature):
        ...
