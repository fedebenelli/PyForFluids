import numpy as np

import pyforfluids.models as models

import pytest


def test_components():
    model = models.GERG2008()
    valid_components = model.valid_components
    model.validate_components(valid_components)

    with pytest.raises(Exception):
        wrong_component = "other_string"
        model.validate_components(valid_components + [wrong_component])


def test_validate_pt_values_min_extended_range():
    random_extended_press = np.random.uniform(35.1, 70.0)
    min_random_extended_temp = np.random.uniform(60.0, 89.9)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(random_extended_press, min_random_extended_temp)


def test_validate_pt_values_max_extended_range():
    random_extended_press = np.random.uniform(35.1, 70.0)
    max_random_extended_temp = np.random.uniform(450.1, 700.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(random_extended_press, max_random_extended_temp)


def test_validate_pt_values_min_invalid_range():
    random_invalid_press = np.random.uniform(70.1, 100.0)
    min_random_invalid_temp = np.random.uniform(00.0, 59.9)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(random_invalid_press, min_random_invalid_temp)


def test_validate_pt_values_max_invalid_range():
    random_invalid_press = np.random.uniform(70.1, 100.0)
    max_random_invalid_temp = np.random.uniform(700.1, 1000.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(random_invalid_press, max_random_invalid_temp)


def test_validate_pt_values_negative_range():
    random_nagetive_press = -np.random.uniform(0.1, 100.0)
    random_negative_temp = -np.random.uniform(0.1, 100.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(random_nagetive_press, random_negative_temp)
