import pyforfluids as pff
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

model = pff.models.GERG2008()
composition = {"ethane": 0.80, "propane": 0.20}

density = 1
pressure = 1.713e6
temperature = 322.04

fluid = pff.Fluid(model, composition, temperature=temperature, density=density)
fluid.set_pressure(pressure)


#x_range = np.linspace(0.01, 0.99, 5)
#
#for x1 in x_range:
#    fluid.composition = {
#            "ethane": x1,
#            "propane": 1-x1
#            }
#    fluid.set_pressure(pressure)
#    res = pff.equilibrium.flash.bub_t(fluid, pressure=pressure, iterations=50)

T, P = pff.equilibrium.flash.envelope(fluid)
