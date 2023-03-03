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
    def ge(cls, z, volume, temperature):
        raise NotImplementedError
