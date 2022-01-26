"""Flash module."""

import numpy as np
from scipy.optimize import root_scalar


def rachford_rice(vapor_fraction, z, K):
    """Rachford rice equation.

    Parameters
    ----------
    vapor_fraction: float
        Vapor fraction.
    z: array
        array of global molar fractions.
    K: array
        array of K-factors.

    Returns
    -------
    g: float
        Valuation of the Rachford rice equation.
    """
    msk = np.where(z != 0)
    g = np.zeros_like(z)
    g[msk] = z[msk] * (K[msk] - 1) / (1 - vapor_fraction + vapor_fraction*K[msk])
    g = g.sum()  

    return g


def update_density(vapor, liquid, pressure):
    rho_l = liquid.density_iterator(pressure)[0]
    rho_g = vapor.density_iterator(pressure)[1]

    liquid.set_density(rho_l)
    vapor.set_density(rho_g)

    liquid.calculate_properties()
    vapor.calculate_properties()

    return vapor, liquid


def update_concentration(x):
    gerg_components = [
        "methane", "nitrogen", "carbon_dioxide",
        "ethane", "propane", "butane",
        "isobutane", "pentane", "isopentane",
        "hexane", "heptane", "octane",
        "nonane", "decane", "hydrogen",
        "oxygen", "carbon_monoxide", "water",
        "hydrogen_sulfide", "helium", "argon"
        ]
    components = {}

    for i, comp in enumerate(gerg_components):
        if x[i] != 0:
            components[comp] = round(x[i], 8)

    return components


def get_vl(fluid, Pesp):
    """
    """
    rho_l, rho_v, _ = fluid.density_iterator(Pesp)

    vapor = fluid.copy()
    vapor.set_density(rho_v)
    vapor.calculate_properties()

    liquid = fluid.copy()
    liquid.set_density(rho_l)
    liquid.calculate_properties()

    return vapor, liquid


def fix_k(z, K, it=0):
    """Increase or decrease the value of K to force the root of the
    Rachford-Rice equation to be in the interval (0, 1).

    Parameters
    ----------
    z: array
        Array of global molar fractions.
    K: array
        Array of K-values.
    it: int
    Number of iteration.

    Returns
    -------
    K: array
        Fixed K value.
    """

    g0 = rachford_rice(1e-9, z, K)
    g1 = rachford_rice(1-1e-9, z, K)
    
    it = it + 1

    if it > 3:
        return K

    if g0 < 0 or g1 > 0:
        while g0 < 0:
            K = 1.1*K
            g0 = rachford_rice(0, z, K)
        while g1 > 0:
            K = 0.9*K
            g1 = rachford_rice(1, z, K)
        return fix_k(z, K, it)
    else:
        return K


def solve_rr(z, K):
    """Solve the Rachford-Rice equation.
    Parameters
    ----------
    z: array
        array of molar fractions
    K: array
        array of K-values

    Returns
    -------
    x: array
        Array of liquid's molar fractions.
    y: array
        Array of vapor's molar fractions.
    vapor_fraction: float
        Vapor fraction.
    """

    b_max = (1-z)/(1-K)
    try:
        b_max = b_max[np.where((b_max <= 1) & (b_max > 0))].min()
    except ValueError:
        b_max = 1

    b_min = (K*z-1)/(K-1)
    try:
        b_min = b_min[np.where((b_min < 1) & (b_min >= 0))].max()
    except ValueError:
        b_min = 0

    # In case there still no solution, assume there the system is single-phase
    g0 = rachford_rice(b_min, z, K)
    g1 = rachford_rice(b_max, z, K)

    if g0 < 0:
        return z, z, 0, 1
    elif g1 > 0:
        return z, z, 1, 1

    root = root_scalar(
            rachford_rice, args=(z, K),
            bracket=(b_min, b_max),
            method='bisect',
            x0=b_min, x1=b_max,
            xtol=1e-4, rtol=1e-4,
            )

    it = root.iterations
    vapor_fraction = root.root

    x = z/(1 - vapor_fraction + vapor_fraction*K)
    y = K*x

    return x, y, vapor_fraction, it


def flash_pt(fluid, pressure, temperature, iterations=50):
    """Vapor-Liquid PT flash.
   
    Parameters
    ----------
    fluid: Fluid
        pyforfluids.Fluid object where to realize calculations.
    pressure: float
        Pressure [Pa].
    temperature: float
        Temperature [K].
    iterations: int
        Max amount of iterations to run.
    """
    z = fluid.model.set_concentration(fluid.composition)
    vapor, liquid = get_vl(fluid, pressure)

    for it in range(iterations):
        k_i = liquid["fugacity_coefficent"] - vapor["fugacity_coefficent"]
        k_i = np.exp(k_i)

        x, y, beta, it_rr = solve_rr(z, k_i)

        #  Update concentrations and then densities
        vapor.composition = update_concentration(y)
        liquid.composition = update_concentration(x)

        # Update the K-values
        k_new = liquid['fugacity_coefficent'] - vapor['fugacity_coefficent']
        k_new = np.exp(k_new)

        if abs((k_new - k_i)/k_i).sum() > 0.2:
            k_new = fix_k(z, k_new)

        if not np.allclose(k_i, k_new, rtol=1e-6, atol=1e-6):
            vapor, liquid = update_density(vapor, liquid, pressure)
            return vapor, liquid, beta, it

    vapor, liquid = update_density(vapor, liquid, pressure)
    return vapor, liquid, beta, it
