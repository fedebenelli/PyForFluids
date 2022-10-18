"""Flash module."""

import numpy as np

from scipy.optimize import root_scalar


def k_wilson(fluid):
    """Wilson K-factor estimation.

    From a given Fluid, get it's K-factors at it's temperature and pressure
    conditions.

    Parameters
    ----------
    fluid: pyforfluids.Fluid
        Fluid to which calculate it's K-factors

    Returns
    -------
    k_i: array_like
        K-factors array
    """
    acentric_factor = fluid.model.w
    critical_pressure = fluid.model.pc
    critical_temperature = fluid.model.tc
    z = fluid.model.set_concentration(fluid.composition)
    pressure = fluid.pressure
    temperature = fluid.temperature

    k_i = np.zeros_like(z)
    msk = np.where(z != 0)

    p_c = critical_pressure[msk]
    t_c = critical_temperature[msk]
    w = acentric_factor[msk]

    k_i[msk] = np.log(p_c / pressure) + 5.373 * (1 + w) * (
        1 - t_c / temperature
    )
    k_i = np.exp(k_i)

    return k_i


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
    """Update two fluids to new compositions and pressure.

    Parameters
    ----------
    liquid: pyforfluids.Fluid
        Liquid fluid
    vapor: pyforfluids.Fluid
        Vapor fluid
    x: array_like
        Array of liquid compositions
    y: array_like
        Array of vapor compositions
    pressure: float
        Pressure

    Returns
    -------
    liquid: pyforfluids.Fluid
        Liquid fluid
    vapor: pyforfluids.Fluid
        Vapor fluid
    """
    liquid = update_concentration(liquid, x)
    vapor = update_concentration(vapor, y)

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


def update_concentration(fluid, x):
    """Update Fluid concentration."""
    components = fluid.model.names
    composition = {}

    for i, comp in enumerate(components):
        if x[i] != 0:
            composition[comp] = round(x[i], 8)

    return fluid.copy(composition=composition)


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

    z = fluid.model.set_concentration(fluid.composition)
    vapor, liquid = get_vl(fluid, pressure)

    k_i = k_wilson(vapor)
    k_i = fix_k(z, k_i)
    x, y, vapor_fraction, it_rr = solve_rr(z, k_i)
    liquid, vapor = update_fluids(liquid, vapor, x, y, pressure)

    for it in range(iterations):
        k_i = np.exp(liquid["lnfug"] - vapor["lnfug"])

        x, y, vapor_fraction, it_rr = solve_rr(z, k_i)

        #  Update concentrations and then densities
        liquid, vapor = update_fluids(liquid, vapor, x, y, pressure)

        # Update the K-values
        k_new = liquid["lnfug"] - vapor["lnfug"]
        k_new = np.exp(k_new)

        if np.allclose(k_i, k_new, rtol=rtol, atol=atol):
            return vapor, liquid, vapor_fraction, it

    liquid, vapor = update_density(liquid, vapor, pressure)
    return vapor, liquid, vapor_fraction, it


# def bub_p(fluid, temperature, iterations=50, rtol=1e-3, atol=1e-3):
#     def p_wilson(temperature, z):
#         w_i = fluid.model.w
#         p_c = fluid.model.pc
#         t_c = fluid.model.tc
#
#         p = z * p_c * np.exp(5.373 * (1 + w_i) * (1 - t_c / temperature))
#         p = np.sum(p)
#
#         return p
#
#     # Make each phase fluid
#     fluid_i = fluid.copy()
#     vapor = fluid.copy()
#     liquid = fluid.copy()
#
#     # Wilson initialization
#     z = fluid.model.set_concentration(fluid.composition)
#     msk = np.where(z != 0)
#     p_i = p_wilson(temperature, z)
#     fluid_i.set_pressure(p_i)
#     k_i = k_wilson(fluid_i)
#     x = z
#     y_i = k_i * x
#     liquid, vapor = update_fluids(liquid, vapor, x, y_i, p_i)
#
#     k_i = np.exp(liquid["lnfug"] - vapor["lnfug"])
#
#     for it in range(1, iterations):
#         liquid, vapor = update_fluids(liquid, vapor, x, y_i, p_i)
#
#         # Define iteration step
#         f = np.dot(k_i, z) - 1
#         dfdp = np.dot(
#             z[msk] * k_i[msk],
#             liquid["dlnfug_dp"][msk] - vapor["dlnfug_dp"][msk],
#         )
#
#         # Iteration step definition
#         step = 1  # it*np.log10(p_i)
#
#         # Iteration is taking too long
#         if 5 < it < 30:
#             step = step * 10
#
#         # NaN pressure reached
#         if np.isnan(p_i - step * f / dfdp):
#             return fluid, fluid, p_i, it
#         p_n = p_i - step * f / dfdp
#
#         # Update fluids
#         liquid, vapor = update_fluids(liquid, vapor, x, y_i, p_n)
#
#         # Define new K and molar fractions
#         k_n = np.exp(liquid["lnfug"] - vapor["lnfug"])
#         x = z
#         y_n = k_n * z
#
#         liquid, vapor = update_fluids(liquid, vapor, x, y_n, p_n)
#
#         # If algo converges, return values, else define new points
#         if np.allclose(p_i, p_n, rtol=rtol, atol=atol) and np.allclose(
#             y_i, y_n, rtol=rtol, atol=atol
#         ):
#             return vapor, liquid, p_n, it
#         else:
#             p_i = p_n
#             y_i = y_n
#
#     return vapor, liquid, p_i, it
#
#
# def bub_t(fluid, pressure, iterations=50, rtol=1e-3, atol=1e-3):
#     """Calculate the bubble temperature and vapor-phase composition.
#
#     Parameters
#     ----------
#     fluid: pyforfluids.core.Fluid
#         Fluid to which calculate it's bubble point.
#     pressure: float
#         Pressure where to calculat the bubble point.
#     iterations: int
#         Max number of iterations to realize.
#     rtol: float
#         Relative tolerance for convergence.
#     atol: float
#         Absolute tolerance for convergence.
#
#     Returns
#     -------
#     vapor: pyforfluids.core.Fluid
#         Vapor phase fluid at bubble conditions.
#     liquid:
#         Liquid phase fluid at bubble conditions.
#     t_n: float
#         Bubble temperature.
#     it: int
#         Number of iterations for convergence.
#     """
#
#     def solve_wilson_temp(z, pressure):
#         """Solve the Wilson's K-factors expression for a given pressure."""
#
#         def f(t, z, pressure):
#             k = k_wilson(fluid)
#             return np.dot(z, k) - 1
#
#         temperature = root_scalar(
#             f, args=(z, pressure), method="brentq", x0=250, bracket=(1, 1000)
#         )
#         return temperature.root
#
#     def new_fluids(liquid, vapor, x, y, pressure, temperature):
#         """Update the fluids"""
#         liquid = update_concentration(liquid, x)
#         vapor = update_concentration(vapor, y)
#
#         liquid.set_temperature(temperature)
#
#         rho_l = liquid.density_iterator(pressure, vapor_phase=False)[0]
#
#         liquid.set_density(rho_l)
#
#         liquid.calculate_properties()
#         vapor.calculate_properties()
#
#         return liquid, vapor
#
#     z = fluid.model.set_concentration(fluid.composition)
#     liquid = fluid.copy()
#     vapor = fluid.copy()
#
#     # Wilson init
#     t_i = solve_wilson_temp(z, pressure)
#     k_i = k_wilson(z, pressure, t_i)
#     x = z
#     y_i = k_i * z
#
#     liquid, vapor = new_fluids(liquid, vapor, x, y_i, pressure, t_i)
#
#     for it in range(1, iterations):
#         f = np.dot(k_i, z) - 1
#         dfdt = (k_i * z * (liquid["dlnfug_dt"] - vapor["dlnfug_dt"])).sum()
#
#         step = 25 * it / t_i
#
#         t_n = t_i - step * f / dfdt
#
#         liquid, vapor = new_fluids(liquid, vapor, x, y_i, pressure, t_n)
#         k_n = np.exp(
#             liquid["lnfug"] - vapor["lnfug"]
#         )
#         y_n = k_n * z
#
#         liquid, vapor = new_fluids(liquid, vapor, x, y_n, pressure, t_n)
#
#         if np.allclose(t_i, t_n, atol=atol) and np.allclose(
#             y_i, y_n, atol=atol
#         ):
#             return vapor, liquid, t_n, it
#         else:
#             k_i = k_n
#             y_i = y_n
#             t_i = t_n
#     return vapor, liquid, t_n, it


# def envelope(fluid):
#     def new_fluids(liquid, vapor, x, y, pressure, temperature):
#         """Update the fluids"""
#         liquid = update_concentration(liquid, x)
#         vapor = update_concentration(vapor, y)
#
#         liquid.set_temperature(temperature)
#
#         rho_l = liquid.density_iterator(pressure, vapor_phase=False)[0]
#
#         liquid.set_density(rho_l)
#
#         liquid.calculate_properties()
#         vapor.calculate_properties()
#
#         return liquid, vapor
#
#     def F(X, liquid, vapor, N, S):
#         """Function to solve"""
#         ns, s = S
#         F = np.zeros(N+2)
#         k = np.exp(X[:N])
#         t = np.exp(X[-2])
#         p = np.exp(X[-1])
#
#         x = Z / (1 - BETA + BETA * k)
#         y = k * x
#
#         liquid, vapor = new_fluids(liquid, vapor, x, y, p, t)
#         fug_diff = (liquid["lnfug"] - vapor["lnfug"])
#
#         for i in range(N):
#             F[i] = X[i] - fug_diff[i]
#         F[N] = (x - y).sum()
#         F[N+1] = X[ns] - s
#         return F, liquid, vapor
#
#     def jacobian(liquid, vapor, temperature, pressure, n_spec, N):
#         Jij = np.zeros((N + 2, N + 2))
#         phi_t_l = liquid["dlnfug_dt"]
#         phi_p_l = liquid["dlnfug_dp"]
#         phi_t_v = vapor["dlnfug_dt"]
#         phi_p_v = vapor["dlnfug_dp"]
#
#         phi_ij_l = liquid["dlnfug_dn"]
#         phi_ij_v = vapor["dlnfug_dn"]
#
#         for i in range(N):
#             for j in range(N):
#                 Jij[i, j] = (
#                     k[j] * Z[j] / (1 - BETA + BETA * k[j]) ** 2
#                     * ((1 - BETA) * phi_ij_v[i][j] + BETA * phi_ij_l[i][j])
#                 )
#                 if i == j:
#                     Jij[i, j] += 1
#                 Jij[N, j] = k[j] * Z[j] / (1 - BETA + BETA * k[j]) ** 2
#
#             Jij[i, N] = temperature * (phi_t_l[i] - phi_t_v[i])
#             Jij[i, N + 1] = pressure * (phi_p_l[i] - phi_p_v[i])
#
#             # TODO: In a future this should be a function
#             Jij[N + 1, n_spec] = 1
#         return Jij
#
#     Z = fluid.model.set_concentration(fluid.composition)
#     BETA = 0
#     N = len(Z)
#     p_i = 0.1e5
#
#     vapor, liquid, t_i, it = bub_t(fluid, p_i)
#
#     # Specification index: N For temperature and N+1 for pressure
#     s = np.log(p_i)
#     n_spec = N + 1
#     delta_spec = np.log(1e3)
#
#     P = []
#     T = []
#
#     for i in range(50):
#         ln_k = (liquid["lnfug"] - vapor["lnfug"])
#         k = np.exp(ln_k)
#         x = np.array([*ln_k, np.log(t_i), np.log(p_i)])
#
#         dx = np.ones_like(x)
#         S = (n_spec, s)
#
#         # Solve F(X) = 0
#         init = 0
#         while abs(dx).max() > 1e-3:
#             init += 1
#             print(i, abs(dx).max())
#             t_i = np.exp(x[-2])
#             p_i = np.exp(x[-1])
#             f, liquid, vapor = F(x, liquid, vapor, N, S)
#             j = jacobian(liquid, vapor, t_i, p_i, n_spec, N)
#
#             dx = np.linalg.solve(j, -f)
#
#             while abs(dx).max() > 5:
#                 dx = dx/100
#
#             if init > 50:
#                 dx *= 100
#             x += dx
#
#         # Jij dXdS + dFdS = 0
#         dXdS = np.linalg.solve(j, j[-1])
#
#         n_old = n_spec
#         n_spec = np.where(dXdS == np.abs(dXdS.max()))[0][0]
#
#         dXdS = dXdS / np.abs(dXdS).max()
#         delta_spec = np.min([1, delta_spec * x[n_spec]/x[n_old]])
#
#         x = x + delta_spec*dXdS
#
#         t_i = np.exp(x[-2])
#         p_i = np.exp(x[-1])
#         k = np.exp(x[:-2])
#
#         S = (n_spec, x[n_spec])
#         P.append(p_i)
#         T.append(t_i)
#         print(P)
#         print(T)
#
#         liquid, vapor = new_fluids(liquid, vapor, Z, k*Z, p_i, t_i)
#     return T, P
