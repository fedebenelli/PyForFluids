"""Cubic EoS models.
"""

from pyforfluids.fortran.fcubiceos import cubic as feos
import pyforfluids.fortran.fcubiceos as fcub
import numpy as np


class CubicEOS:
    """Cubic EoS object.

    Represents a cubic equation of state model.

    Attributes
    ----------
    model: str
        Equation of State to be used. Available models:
            - 'PR': PengRobinson Equation of State
            - 'SRK': SRK Equation of State
            - 'RKPR': RKPR Equation of State
    critical_atractive: array_like
        Critical attractive parameter (ac)
    repulsive_parameter: array_like
        Repulsive parameter (b)
    delta_1_parameter: array_like
        delta_1 parameter, defaults to 1 + sqrt(2) for PengRobinson EoS and 1
        for SRK EoS. Must be specified for the RKPR EoS
    critical_k_parameter: array_like
        K constant used for the calculation of the attractive parameter as
        a function of T.
    critical_temperature: array_like
        Array of critical temperatures.
    critical_pressure: array_like
        Array of critical pressures.
    critical_density: array_like
        Array of critical densities.
    accentric_factor: array_like
        Array of accentric factors.
    mix_rule: str
        Mixing rule of fluids.
            - 'quadratic': classic quadratic mixing rules.
    kij0_matrix: 2D-array_like
        Constant binary interaction parameters.
    kijinf_matrix: 2D-array_like
        Kij_inf matrix that's used to calculate the binary interaction
        parameters as a function of temperature.
    t_star_matrix
        T_star matrix that's used to calculate the binary interaction
        parameters as a function of temperature.
    lij_matrix: 2D-array_like
        Binary repulsion parameters.
    volume_traslation: array_like
    volume_shift: array_like

    Methods
    -------
    validate_components:
        For Fluid compatibility, always returns True.
    set_concentration:
        Get a concentrations vector from the composition.
    calculate_properties:
        Calculate the available properties.
    """

    # Gas constant from inside the Fortran code
    R = feos.rgas

    name = "Cubic EoS"

    # Selectors to setup the Fortran commons values later
    __model_selector = {"SRK": 1, "PR": 2, "RKPR": 3}
    __mixing_rules_selector = {"quadatric": 0, "cubic": 3}

    def __init__(
        self,
        model,
        critical_atractive,
        repulsive_parameter,
        delta_1_parameter,
        critical_k_parameter,
        critical_temperature,
        critical_pressure,
        critical_density,
        accentric_factor,
        mix_rule,
        kij0_matrix,
        kijinf_matrix,
        t_star_matrix,
        lij_matrix,
        volume_traslation,
        volume_shift,
    ):
        self.model = model
        self.nmodel = self.__model_selector[model]
        self.ncomb = self.__mixing_rules_selector[mix_rule]

        self.ac = critical_atractive
        self.b = repulsive_parameter
        self.k = critical_k_parameter

        self.tc = critical_temperature
        self.pc = critical_pressure
        self.w = accentric_factor

        self.kij = kij0_matrix
        self.lij = lij_matrix

        if model == "RKPR":
            self.delta_1 = delta_1_parameter
        elif model == "PR":
            self.delta_1 = 1 + np.sqrt(2)
        elif model == "SRK":
            self.delta_1 = 1

        if kijinf_matrix is not None:
            self.tdep = 1
            self.kijinf_matrix = kijinf_matrix
            self.t_star_matrix = t_star_matrix
        else:
            self.tdep = 0
            self.kijinf_matrix = 0
            self.t_star_matrix = 0

        if volume_traslation is not None:
            self.volume_traslation = volume_traslation
            self.volume_shift = volume_shift
        else:
            self.volume_traslation = 0
            self.volume_shift = 0

    def validate_components(self, composition):
        return True

    def set_concentration(self, composition):
        return np.array([composition[i] for i in composition], dtype="d")

    def calculate_properties(
        self, temperature, pressure, density, composition
    ):
        """Calculate the thermodynamic properties of the given fluid.

        Calculation of the thermodynamic properties of the fluid at it's given
        temperature and density.

        Parameters
        ----------
        temperature: float
            Fluid temperature in Kelvin degrees [K]
        pressure: float
            Fluid pressure in Pascal [Pa]
        density: float
            Fluid density in mol per liter [mol/L]
        composition: dict
            Dictionary of the compounds concentrations as:
            ``{"methane": 0.8, "ethane":0.2}``
            When necessary, the concentration values are normalized.

        Returns
        -------
        dict
            Dictionary of the thermodynamic properties of the given fluid.
        """
        self.__set_commons()
        R = self.R

        volume = 1 / density
        concentrations = self.set_concentration(composition)
        number_of_moles = concentrations.sum()

        residual_helmholtz_and_derivatives = feos.helmholtz(
            concentrations, volume, temperature
        )

        dar_dv = residual_helmholtz_and_derivatives[1]
        dar_dv2 = residual_helmholtz_and_derivatives[3]
        dar_dn = residual_helmholtz_and_derivatives[4]

        # Properties
        p = R * temperature / volume - dar_dv
        z = (p * volume) / (R * temperature)

        lnfug = dar_dn - np.log(z)

        dp_drho = dar_dv2 * volume**2 + number_of_moles * R * temperature

        return {
            "pressure": p,
            "residual_helmholtz": residual_helmholtz_and_derivatives[0],
            "dar_dv": dar_dv,
            "dar_dtv": residual_helmholtz_and_derivatives[2],
            "compressibility_factor": z,
            "dp_drho": dp_drho,
            "lnfug": lnfug,
        }

    def __set_commons(self):
        """Set up the Fortran model commons."""
        nc = len(self.ac)

        # Model to be used
        fcub.model.nmodel = self.nmodel

        # Components parameters
        fcub.components.ac[:nc] = self.ac
        fcub.components.b[:nc] = self.b
        fcub.components.rm[:nc] = self.k
        fcub.components.ntdep = self.tdep
        fcub.components.kij[:nc, :nc] = self.kij
        fcub.components.del1[:nc] = self.delta_1

        # Critical constants
        fcub.crit.tc[:nc] = self.tc
        fcub.crit.pc[:nc] = self.pc
        fcub.crit.dceos[:nc] = 0  # Variable not used in the original code
        fcub.crit.om[:nc] = self.w

        # Combining rule
        fcub.rule.ncomb = self.ncomb

        # Temperature dependant binary interaction parameters
        fcub.tdep.kinf[:nc, :nc] = self.kijinf_matrix
        fcub.tdep.tstar[:nc, :nc] = self.t_star_matrix

        # Binary repulsive parameters
        fcub.lforin.lij[:nc, :nc] = self.lij

        # Repulsive parameter matrix
        b = np.array([self.b])
        bij = (1 - self.lij) * (b + b.T) / 2
        fcub.bcross.bij[:nc, :nc] = bij

    def __check_model(self, model):
        """Check if the model is available.

        Parameters
        ----------
        model: str
            Model to be used, can be checked with:
                `pyforfluids.models.CubicEOS.model_selector`

        Returns
        -------
        bool
            Is the selected model available?
        """
        return model in self.__model_selector

    def __repr__(self):
        """Model representation."""
        return f"{self.name}: {self.model}"
