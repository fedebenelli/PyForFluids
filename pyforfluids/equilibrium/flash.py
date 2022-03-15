"""Flash module."""

import numpy as np

from scipy.optimize import root_scalar

from ..fortran import gerg2008f
from ..core import Fluid


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


def update_fluids(liquid, vapor, x, y, pressure):
    liquid.composition = update_concentration(x)
    vapor.composition = update_concentration(y)
    liquid, vapor = update_density(liquid, vapor, pressure)

    liquid.calculate_properties()
    vapor.calculate_properties()

    return liquid, vapor


def update_density(liquid, vapor, pressure):
    """Update both phases densities."""
    rho_l = liquid.density_iterator(pressure, vapor_phase=False)[0]
    rho_g = vapor.density_iterator(pressure, liquid_phase=False)[1]

    liquid.set_density(rho_l)
    vapor.set_density(rho_g)

    return liquid, vapor


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
    fluid, pressure, temperature, iterations=50, rtol=1e-5, atol=1e-5
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

    z = fluid.model.normalize(fluid.composition)
    vapor, liquid = get_vl(fluid, pressure)

    k_i = k_wilson(z, pressure, temperature)
    k_i = fix_k(z, k_i)
    x, y, vapor_fraction, it_rr = solve_rr(z, k_i)
    liquid, vapor = update_fluids(liquid, vapor, x, y, pressure)

    for it in range(iterations):
        k_i = np.exp(
            liquid["lnfug"] - vapor["lnfug"]
        )

        x, y, vapor_fraction, it_rr = solve_rr(z, k_i)

        #  Update concentrations and then densities
        liquid, vapor = update_fluids(liquid, vapor, x, y, pressure)

        # Update the K-values
        k_new = liquid["lnfug"] - vapor["lnfug"]
        k_new = np.exp(k_new)

        #while abs((k_new - k_i) / k_i).sum() > 0.2:
        #    k_new = fix_k(z, k_new)

        if np.allclose(k_i, k_new, rtol=rtol, atol=atol):
            return vapor, liquid, vapor_fraction, it

    liquid, vapor = update_density(liquid, vapor, pressure)
    return vapor, liquid, vapor_fraction, it


def k_wilson(z, pressure, temperature):
    acentric_factor = gerg2008f.parameters.acentric_factor
    critical_pressure = gerg2008f.parameters.p_c
    critical_temperature = gerg2008f.parameters.t_c

    K_i = np.zeros_like(z)
    msk = np.where(z != 0)

    p_c = critical_pressure[msk]
    t_c = critical_temperature[msk]
    w = acentric_factor[msk]

    K_i[msk] = np.log(p_c / pressure) + 5.373 * (1 + w) * (
        1 - t_c / temperature
    )
    K_i = np.exp(K_i)

    return K_i


def p_wilson(temperature, z):
    w_i = gerg2008f.parameters.acentric_factor
    p_c = gerg2008f.parameters.p_c
    t_c = gerg2008f.parameters.t_c

    p = z * p_c * np.exp(5.373 * (1 + w_i) * (1 - t_c / temperature))
    p = np.sum(p)

    return p


def bub_p(fluid, temperature, iterations=50, rtol=1e-5, atol=1e-5):
    # Make each phase fluid
    vapor = fluid.copy()
    liquid = fluid.copy()

    # Wilson initialization
    z = fluid.model.set_concentration(fluid.composition)
    msk = np.where(z != 0)
    p_i = p_wilson(temperature, z)
    k_i = k_wilson(z, p_i, temperature)
    x = z
    y_i = k_i * x

    for it in range(1, iterations):
        liquid, vapor = update_fluids(liquid, vapor, x, y_i, p_i)
        # Define iteration step
        f = np.dot(k_i, z) - 1
        dfdp = np.dot(
                z[msk] * k_i[msk],
                liquid["dlnfug_dp"][msk] - vapor["dlnfug_dp"][msk]
                )

        # Iteration step definition
        step = 1e5 * it

        # This is taking too long
        if it > 5:
            step = step * 10

        p_n = p_i - step * f / dfdp

        # NaN pressure reached, no convergence so assume there is no
        # separation here
        if np.isnan(p_n):
            return fluid, fluid, p_i, it

        # Update fluids
        liquid, vapor = update_fluids(liquid, vapor, x, y_i, p_n)

        # Define new K and molar fractions
        k_n = np.exp(
            liquid["fugacity_coefficent"] - vapor["fugacity_coefficent"]
        )
        x = z
        y_n = k_n * z

        liquid, vapor = update_fluids(liquid, vapor, x, y_n, p_n)
        print(
            liquid.composition["ethane"],
            vapor.composition["ethane"],
            it,
            p_n,
            k_n[3],
        )

        # If algo converges, return values, else define new points
        if np.allclose(y_i, y_n, rtol=rtol, atol=atol) and np.allclose(
            p_i, p_n, rtol=rtol, atol=atol
        ):
            return vapor, liquid, p_n, it
        else:
            p_i = p_n
            y_i = y_n

    return vapor, liquid, p_i, it


def envelope(fluid):
    Z = fluid.model.set_concentration(fluid.composition)
    msk = np.where(Z != 0)
    Z = Z[msk]
    BETA = 0
    vapor, liquid, p_ini, _ = bub_p(fluid, 200, 10, 1e-3, 1e-3)

    # Specification index: N For temperature and N+1 for pressure
    n_spec = N + 1 

    def jac(liquid, vapor, temperature, pressure, n_spec, N):
        Jij = np.zeros_like(N + 2, N + 2)
        phi_t_l = liquid["dlnfug_dt"]
        phi_p_l = liquid["dlnfug_dp"]
        phi_t_v = vapor["dlnfug_dt"]
        phi_p_v = vapor["dlnfug_dp"]

        k = np.exp(
            liquid["fugacity_coefficent"] - vapor["fugacity_coefficent"]
        )

def envelope(fluid):
    def new_fluids(liquid, vapor, x, y, pressure, temperature):
        """Update the fluids"""
        __import__('ipdb').set_trace()
        liquid_composition = update_concentration(x)
        vapor_composition = update_concentration(y)

        liquid.set_composition(liquid_composition)
        liquid.set_temperature(temperature)

        vapor = Fluid(
            model=vapor.model,
            composition=vapor_composition,
            pressure=pressure,
            temperature=temperature,
        )

        rho_l = liquid.density_iterator(pressure, vapor_phase=False)[0]

        liquid.set_density(rho_l)

        liquid.calculate_properties()
        vapor.calculate_properties()

        return liquid, vapor
    
    def F(X, liquid, vapor, N):
        F = np.zeros(N+2)
        k = np.exp(X[:N])
        t = np.exp(X[-2])
        p = np.exp(X[-1])

        x = Z / (1 - BETA + BETA * k)
        y = k * x
        liquid, vapor = new_fluids(liquid, vapor, x, y, p, t)

        fug_diff = (liquid["lnfug"] - vapor["lnfug"])[msk]
        for i in range(N):
            F[i] = X[i] - fug_diff[i]
        F[N] = rachford_rice(BETA, Z, X[:N])
        F[N+1] = 0
        return F

    def jacobian(liquid, vapor, temperature, pressure, n_spec, N):
        Jij = np.zeros((N + 2, N + 2))
        phi_t_l = liquid["dlnfug_dt"][msk]
        phi_p_l = liquid["dlnfug_dp"][msk]
        phi_t_v = vapor["dlnfug_dt"][msk]
        phi_p_v = vapor["dlnfug_dp"][msk]

        phi_ij_l = liquid["dlnfug_dn"][msk]
        phi_ij_v = vapor["dlnfug_dn"][msk]

        for i in range(N):
            for j in range(N):
                Jij[i, j] = (
                    k[j] * Z[j]
                    / (1 - BETA + BETA * k[j]) ** 2
                    * (
                        (1 - BETA) * phi_ij_v[i][msk][j] 
                        + BETA * phi_ij_l[i][msk][j]
                        )
                )
                if i == j:
                    Jij[i, j] += 1
                Jij[N, j] = k[j] * Z[j] / (1 - BETA + BETA * k[j]) ** 2

            Jij[i, N] = temperature * (phi_t_l[i] - phi_t_v[i])
            Jij[i, N + 1] = pressure * (phi_p_l[i] - phi_p_v[i])

            # TODO: In a future this should be a function
            Jij[N + 1, n_spec] = 1
        return Jij

    Z = fluid.model.normalize(fluid.composition)
    msk = np.where(Z != 0)
    #Z = Z[msk]
    BETA = 0
    N = len(Z)
    p_i = 0.1e6

    vapor, liquid, t_i, it = bub_t(fluid, p_i)

    # Specification index: N For temperature and N+1 for pressure
    n_spec = N + 1
    delta_spec = np.log(1e5)

    for i in range(50):
        ln_k = (liquid["lnfug"] - vapor["lnfug"])[msk]
        k = np.exp(ln_k)
        x = np.array([*ln_k, np.log(t_i), np.log(p_i)])

        j = jacobian(liquid, vapor, t_i, p_i, n_spec, N)
        f = F(x, liquid, vapor, N)

        # Jij dXdS + dFdS = 0
        dXdS = np.linalg.solve(j, j[-1])
        n_old = n_spec
        n_spec = np.where(dXdS == np.abs(dXdS.max()))[0][0]
        dXdS = dXdS / np.abs(dXdS).max()
        delta_spec = delta_spec * x[n_spec]/x[n_old]

        x = x + delta_spec*dXdS
