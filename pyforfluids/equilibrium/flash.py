"""Flash module."""

import numpy as np

from scipy.optimize import root_scalar

from ..fortran import gerg2008f


def rachford_rice(vapor_fraction, z, k):
    """Rachford rice equation.

    Parameters
    ----------
    vapor_fraction: float
        Vapor fraction.
    z: array
        array of global molar fractions.
    k: array
        array of k-factors.

    Returns
    -------
    g: float
        Valuation of the Rachford rice equation.
    """
    msk = np.where(z != 0)
    g = np.zeros_like(z)
    g[msk] = (
        z[msk] * (k[msk] - 1) / (1 - vapor_fraction + vapor_fraction * k[msk])
    )
    g = g.sum()

    return g


def update_density(vapor, liquid, pressure):
    """Update both phases densities."""
    rho_l = liquid.density_iterator(pressure, vapor_phase=False)[0]
    rho_g = vapor.density_iterator(pressure, liquid_phase=False)[1]

    liquid.set_density(rho_l)
    vapor.set_density(rho_g)

    liquid.calculate_properties()
    vapor.calculate_properties()

    return vapor, liquid


def update_concentration(x):
    """Update Fluid concentration."""
    gerg_components = [
        "methane",
        "nitrogen",
        "carbon_dioxide",
        "ethane",
        "propane",
        "butane",
        "isobutane",
        "pentane",
        "isopentane",
        "hexane",
        "heptane",
        "octane",
        "nonane",
        "decane",
        "hydrogen",
        "oxygen",
        "carbon_monoxide",
        "water",
        "hydrogen_sulfide",
        "helium",
        "argon",
    ]
    components = {}

    for i, comp in enumerate(gerg_components):
        if x[i] != 0:
            components[comp] = round(x[i], 8)

    return components


# TODO: This method could be part of fluid or a separate module, including
# phase stability analysis.
def get_vl(fluid, pressure):
    """Get vapor and liquid Fluids based on a pressure."""
    rho_l, rho_v, _ = fluid.density_iterator(pressure)

    vapor = fluid.copy()
    vapor.set_density(rho_v)
    vapor.calculate_properties()

    liquid = fluid.copy()
    liquid.set_density(rho_l)
    liquid.calculate_properties()

    return vapor, liquid


def fix_k(z, k, it=0):
    """Fix K-values to assure convergence.

    Increase or decrease the value of k to force the root of the
    Rachford-Rice equation to be in the interval (0, 1).

    Parameters
    ----------
    z: array
        Array of global molar fractions.
    k: array
        Array of k-values.
    it: int
    Number of iteration.

    Returns
    -------
    k: array
        Fixed k value.
    """
    g0 = rachford_rice(1e-9, z, k)
    g1 = rachford_rice(1 - 1e-9, z, k)

    it = it + 1

    if it > 3:
        return k

    if g0 < 0 or g1 > 0:
        while g0 < 0:
            k = 1.1 * k
            g0 = rachford_rice(0, z, k)
        while g1 > 0:
            k = 0.9 * k
            g1 = rachford_rice(1, z, k)
        return fix_k(z, k, it)
    else:
        return k


def solve_rr(z, k):
    """Solve the Rachford-Rice equation.

    Parameters
    ----------
    z: array
        array of molar fractions
    k: array
        array of k-values

    Returns
    -------
    x: array
        Array of liquid's molar fractions.
    y: array
        Array of vapor's molar fractions.
    vapor_fraction: float
        Vapor fraction.
    """
    # Find limits of vapor_fraction where to find the root
    b_max = (1 - z) / (1 - k)
    try:
        b_max = b_max[np.where((b_max <= 1) & (b_max > 0))].min()
    except ValueError:
        b_max = 1

    b_min = (k * z - 1) / (k - 1)
    try:
        b_min = b_min[np.where((b_min < 1) & (b_min >= 0))].max()
    except ValueError:
        b_min = 0

    # In case there still no solution, assume there the system is single-phase
    g0 = rachford_rice(b_min, z, k)
    g1 = rachford_rice(b_max, z, k)

    if g0 < 0:
        return z, z, 0, 1
    elif g1 > 0:
        return z, z, 1, 1

    root = root_scalar(
        rachford_rice,
        args=(z, k),
        bracket=(b_min, b_max),
        method="bisect",
        x0=b_min,
        x1=b_max,
        xtol=1e-4,
        rtol=1e-4,
    )

    it = root.iterations
    vapor_fraction = root.root

    x = z / (1 - vapor_fraction + vapor_fraction * k)
    y = k * x

    return x, y, vapor_fraction, it


def flash_pt(
    fluid, pressure, temperature, iterations=50, rtol=1e-10, atol=1e-10
):
    """Vapor-Liquid PT flash.

    Solve the Flash-PT for a given Fluid at given pressure and temperature.

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
    rtol: float
        Relative tolerance between K-values to accept convergence.
    atol: float
        Absolute tolerance between K-values to accept convergence.

    Returns
    -------
    vapor: Fluid
        Fluid object with the composition of the vapor-phase.
    liquid: Fluid
        Fluid object with the composition of the liquid-phase.
    vapor_fraction: float
        Fraction of vapor in the mixture.
    """
    z = fluid.model.set_concentration(fluid.composition)
    vapor, liquid = get_vl(fluid, pressure)

    for it in range(iterations):
        it += 1
        k_i = liquid["fugacity_coefficent"] - vapor["fugacity_coefficent"]
        k_i = np.exp(k_i)

        x, y, vapor_fraction, it_rr = solve_rr(z, k_i)

        #  Update concentrations and then densities
        vapor.composition = update_concentration(y)
        liquid.composition = update_concentration(x)
        vapor, liquid = update_density(vapor, liquid, pressure)

        # Update the K-values
        k_new = liquid["fugacity_coefficent"] - vapor["fugacity_coefficent"]
        k_new = np.exp(k_new)

        if abs((k_new - k_i) / k_i).sum() > 0.2:
            k_new = fix_k(z, k_new)

        if np.allclose(k_i, k_new, rtol=rtol, atol=atol):
            return vapor, liquid, vapor_fraction, it

    vapor, liquid = update_density(vapor, liquid, pressure)
    return vapor, liquid, vapor_fraction, it


def K_wilson(z, pressure, temperature):
    acentric_factor = gerg2008f.parameters.acentric_factor
    critical_pressure = gerg2008f.parameters.p_c
    critical_temperature = gerg2008f.parameters.t_c

    n = len(z)
    K_i = np.full(n, -np.inf, dtype='d')
    K_i = np.nan_to_num(K_i)
    msk = np.where(z != 0)

    p_c = critical_pressure[msk]
    t_c = critical_temperature[msk]
    w = acentric_factor[msk]

    K_i[msk] = np.log(p_c/pressure) + 5.373*(1+w)*(1-t_c/temperature)
    K_i = np.exp(K_i)

    return K_i


def p_wilson(temperature, z):
    w_i = gerg2008f.parameters.acentric_factor
    p_c = gerg2008f.parameters.p_c
    t_c = gerg2008f.parameters.t_c

    p = z*p_c*np.exp(5.373*(1+w_i)*(1-t_c/temperature))
    p = np.sum(p)

    return p


def bub_p(fluid, temperature, iterations=50, rtol=1e-10, atol=1e-10):
    z = fluid.model.set_concentration(fluid.composition)
    p_i = p_wilson(temperature, z)
    k_i = K_wilson(z, p_i, temperature)

    vapor = fluid.copy()
    liquid = fluid.copy()

    x = z
    y = k_i*x

    liquid.composition = update_concentration(x)
    vapor.composition = update_concentration(y)
    vapor, liquid = update_density(vapor, liquid, p_i)

    k_i = np.exp(liquid['fugacity_coefficent'] - vapor['fugacity_coefficent'])

    for it in range(iterations):
        f = np.dot(k_i, z) - 1
        dfdp = z * k_i * (liquid['dlnfug_dp'] - vapor['dlnfug_dp'])
        dfdp = np.sum(dfdp)

        step = 1
        while abs(step*f/dfdp) < p_i/10000:
            step = step*5

        p_i = p_i - step*f/dfdp
        print(liquid.composition['ethane'], vapor.composition["ethane"], it, p_i)

        liquid.composition = update_concentration(x)
        vapor.composition = update_concentration(y)
        vapor, liquid = update_density(vapor, liquid, p_i)
        
        vapor.calculate_properties()
        liquid.calculate_properties()

        k_new = np.exp(
                liquid['fugacity_coefficent'] - vapor['fugacity_coefficent']
        )

        if np.allclose(k_i, k_new, rtol=rtol, atol=atol):
            return vapor, liquid, p_i, it
        else:
            k_i = k_new

        x = z
        y = k_i*z
    return vapor, liquid, p_i, it
