"""
"""

from ABC import ABCMeta, abstractmethod


class ArModel(ABCMeta):

    @abstractmethod
    def residual_helmholtz(self, z, v, t):
        """Method to calculate residual helmholtz
        """
        raise NotImplementedError
