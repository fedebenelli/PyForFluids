#!/usr/bin/env python
# -*- coding: utf-8 -*-
# License: MIT License (https://tldrlegal.com/license/mit-license)
# Copyright (c) 2021-2022 Federico Benelli and others.
# All rights reserved.
"""GERG2008 EoS."""


class GERG2008:
    """GERG2008 equation of state.

    GERG2008 Equation of state described by O. Kunz and W. Wagner [1]_

    The components must be those of the GERG2008 model.
    This class use imported methods from Fortran subroutines for high speed
    calculation of properties.

    Methods
    -------

    References
    ----------
    .. [1] O. Kunz and W. Wagner,
       "The GERG-2008 Wide-Range Equation of State for Natural Gases and
       Other Mixtures: An Expansion of GERG-2004", J. Chem. Eng. Data 2012,
       57, 11, 3032–3091. doi:10.1021/je300655b.
       `<https://pubs.acs.org/doi/10.1021/je300655b>`_
    """
