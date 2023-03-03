"""
"""

from ABC import ABCMeta, abstractmethod


class ArModel(ABCMeta):
    """
    """
    @abstractmethod
    def residual_helmholtz(self, z, volume, temperature):
        """Method to calculate residual helmholtz
        """
        raise NotImplementedError


class GeModel(ABCMeta):
    """
    """
    @abstractmethod
    def ge(self, z, volume, temperature):
        raise NotImplementedError

