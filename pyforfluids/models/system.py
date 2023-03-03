"""
"""
from jax import jacrev, jit


class System:
    def __init__(self, armodel, gemodel, idealmodel):
        self.armodel = armodel
        self.gemodel = gemodel
        self.idealmodel = idealmodel

        if armodel.autodiff:
            self._ar = jit(armodel.residual_helmholtz)
            self._dar = jit(
                jacrev(armodel.residual_helmholtz, argnums=[0, 1, 2])
            )

    def residual_helmholtz(self, z, volume, temperature):
        ar = self._ar(z, volume, temperature)
        dar = self._dar(z, volume, temperature)
        
        return ar, dar

    def ln_phi(self, z, volume, temperature):
        # residual_helmholtz = self.residual_helmholtz(z, volume, temperature)
        ...

    def pressure(self, z, volume, temperature):
        ...

    def volume(self, z, pressure, temperature):
        ...