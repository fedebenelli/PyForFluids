#!/usr/bin/env python
# -*- coding: utf-8 -*-
# License: MIT License (https://tldrlegal.com/license/mit-license)
# Copyright (c) 2021-2022 Federico Benelli and others.
# All rights reserved.

"""PyForFluids.

Fluid properties simulation based on Ecuations of State.
"""

__version__ = "0.0.1a3"

from .core import Fluid  # noqa
from . import fortran  # noqa
from . import models  # noqa
