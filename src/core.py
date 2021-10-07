"""PyForFluids
   ===========

   Module that does stuff.....
"""


class Fluid:
    def __init__(
        self, model, composition, temperature, pressure=None, density=None
    ):
        if (pressure, density) == (None, None) or None not in (
            pressure,
            density,
        ):
            raise ValueError("Presión y densidad no pueden ser ambos None....")

        self.composition = composition
        self.model = model
        self.model.validate_components(self.composition)
        self.temperature = temperature
        self.pressure = pressure
        self.density = density

    def set_temperature(self, temperature):
        self.temperature = temperature

    def set_density(self, density):
        self.density = density

    def set_pressure(self, pressure):
        self.pressure = pressure

    def calculate_properties(self):
        self.properties = self.model.calculate_properties(
            self.temperature, self.pressure, self.density, self.composition
        )
