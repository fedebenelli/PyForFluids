"""
"""
from jax import jit, jacrev
import jax.numpy as np


class System:

    def __init__(self, armodel, gemodel, idealmodel):
        self.armodel = armodel
        self.gemodel = gemodel
        self.idealmodel = idealmodel

        if armodel.autodiff:
            self._residual_helmholtz = jit(
                    jacrev(armodel.residual_helmholtz)
                )

    def residual_helmholtz(self, z, v, t):
        return self._residual_helmholtz(z, v, t)

    def ln_phi(self, z, volume, temperature):
        residual_helmholtz = self.residual_helmholtz(z, volume, temperature)
        ...        
    
    def pressure(self, z, volume, temperature):
        ...

    def volume(self, z, pressure, temperature):
        ...
