"""Cubic EoS models
"""

from pyforfluids.fortran.fcubiceos import cubic as feos
import numpy as np


class CubicEOS:
    """Cubic EoS
    """

    model_selector = {"SRK": 1, "PR": 2, "RKPR": 3}
    mixing_rules_selector = {"quadatric": 0, "cubic": 3}

    def __init__(
            self, model,
            critical_atractive, repulsive_parameter, delta_1_parameter, 
            critical_k_parameter, critical_temperature, critical_pressure, 
            critical_density, accentric_factor, mix_rule, kij0_matrix,
            kijinf_matrix, t_star_matrix, lij_matrix, 
            volume_traslation, volume_shift
            ):
        self.nmodel = self.model_selector[model]
        self.ncomb = self.mixing_rules_selector[mix_rule]

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

    def _setup(self):
        """Set up the Fortran model.
        """
        nc = len(self.ac)
        fcub.model.nmodel = self.nmodel

        fcub.components.ac[:nc] = self.ac
        fcub.components.b[:nc] = self.b
        fcub.components.rm[:nc] = self.k
        fcub.components.ntdep = self.tdep
        fcub.components.kij[:nc, :nc] = self.kij
        fcub.components.del1[:nc] = self.delta_1
        fcub.crit.tc[:nc] = self.tc
        fcub.crit.pc[:nc] = self.pc
        fcub.crit.dceos[:nc] = 0
        fcub.crit.om[:nc] = self.w
        fcub.rule.ncomb = self.ncomb
        fcub.tdep.kinf[:nc, :nc] = 0
        fcub.tdep.tstar[:nc, :nc] = 0
        fcub.lforin.lij[:nc, :nc] = self.lij

        bij = np.zeros((nc, nc))

        for i in range(nc):
            for j in range(i, nc):
                bij[i, j] = (1 - self.lij[i, j]) * (self.b[i] + self.b[j])/2

        fcub.bcross.bij[:nc, :nc] = bij

    def check_model(self, model):
        """Check if the model is available.

        Parameters
        ----------

        model: str
            Model to be used, can be checked with:
                `pyforfluids.models.CubicEOS.model_selector`

        Returns
        -------

        is_available: bool
            Is the selected model available?
        """
        return model in self.model_selector

    def set_concentration(self, composition):
        return np.array(
                [composition[i] for i in composition]
        )

    def calculate_properties(
        self, temperature, pressure, density, composition
    ):
        volume = 1/density
        indicator = 4

        concentrations = self.set_concentration(composition)

        pressure = feos.pressure_calc(concentrations, volume, temperature)

        liquid_root = feos.lnfug(
                1, indicator, temperature, pressure, concentrations
        )
        vapor_root = feos.lnfug(
                -1, indicator, temperature, pressure, concentrations
        )
        lower_g_root = feos.lnfug(
                0, indicator, temperature, pressure, concentrations
        )

        return {
                "pressure": pressure,
                "liquid_root": liquid_root,
                "vapor_root": vapor_root,
                "lower_g_root": lower_g_root
        }
