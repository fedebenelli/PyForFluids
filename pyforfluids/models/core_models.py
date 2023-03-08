"""
Abstract objects used to define the basic skeletons needed.
"""

from abc import ABCMeta, abstractmethod


class ArModel(metaclass=ABCMeta):
    """Residual Helmholtz Model"""

    @abstractmethod
    def residual_helmholtz(self, z, volume, temperature):
        """Method to calculate residual helmholtz"""
        raise NotImplementedError


class GeModel(metaclass=ABCMeta):
    """Excess Gibbs Model"""

    @abstractmethod
    def excess_gibbs(self, z, volume, temperature):
        """Method to calculate excess Gibbs energy"""
        raise NotImplementedError


class CubicMixingRule(metaclass=ABCMeta):
    """Mixing Rule to use inside a Cubic Equation of State."""

    @abstractmethod
    def mix_a(self, z, v, t, a):
        raise NotImplementedError

    @abstractmethod
    def mix_b(self, z, v, t, b):
        raise NotImplementedError

    @abstractmethod
    def mix_c(self, z, v, t, c):
        raise NotImplementedError

    @abstractmethod
    def mix_delta1(self, z, v, t, delta1):
        raise NotImplementedError

    @abstractmethod
    def mix_delta2(self, z, v, t, delta2):
        raise NotImplementedError
