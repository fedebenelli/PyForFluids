import warnings


class InvalidWorkingConditionsValueError(ValueError):
    def __init__(
        self, message="Working conditions belong to the invalid range."
    ):
        self.message = message
        super().__init__(self.message)


def pt_condition_validator(pressure, temperature, show_plot=False):

    if (
        temperature > 700.0
        or temperature < 60.0
        or pressure > 70.0
        or pressure <= 0.0
    ):
        raise InvalidWorkingConditionsValueError()
    elif (
        temperature <= 90.0
        or temperature <= 450.0
        and pressure <= 0.0
        or pressure <= 35.0
    ):
        return print(
            "Working conditions belong to the normal range of validity."
        )
    else:
        warnings.warn(
            "Working conditions belong to the extended vality range."
        )
        return
