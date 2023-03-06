"""
"""

from abc import ABCMeta, abstractmethod


class ArModel(ABCMeta):
    """ """

    @abstractmethod
    def residual_helmholtz(cls, z, volume, temperature):
        """Method to calculate residual helmholtz"""
        raise NotImplementedError


class GeModel(ABCMeta):
    """ """

    @abstractmethod
    def excess_gibbs(cls, z, volume, temperature):
        """Method to calculate excess Gibbs energy"""
        raise NotImplementedError
