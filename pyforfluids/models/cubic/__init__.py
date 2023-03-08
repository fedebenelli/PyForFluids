"""Cubic Equations of State

This module contains all the implemented Cubic Equations of State and also the
respective base class that can be used as an skeleton to implement new models
"""

__all__ = ["CubicEoS", "PR76"]

from pyforfluids.models.cubic.cubic_eos import CubicEoS
from pyforfluids.models.cubic.peng_robinson_76 import PR76
