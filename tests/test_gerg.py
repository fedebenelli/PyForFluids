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


# -- Valid ranges tests -------------------------------------------------------


@pytest.mark.skip(reason="no way of currently testing this")
def test_validate_pt_values_no_warnings():
    random_normal_press = np.random.uniform(0.0, 35.0)
    random_normal_temp = np.random.uniform(90.0, 450.0)
    model = models.GERG2008()
    model.validate_ranges(random_normal_temp, random_normal_press)


def test_validate_pt_values_min_extended_range():
    random_extended_press = np.random.uniform(35.1, 70.0)
    min_random_extended_temp = np.random.uniform(60.0, 89.9)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(min_random_extended_temp, random_extended_press)


def test_validate_pt_values_max_extended_range():
    random_extended_press = np.random.uniform(35.1, 70.0)
    max_random_extended_temp = np.random.uniform(450.1, 700.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(max_random_extended_temp, random_extended_press)


def test_validate_pt_values_min_invalid_range():
    random_invalid_press = np.random.uniform(70.1, 100.0)
    min_random_invalid_temp = np.random.uniform(00.0, 59.9)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(min_random_invalid_temp, random_invalid_press)


def test_validate_pt_values_max_invalid_range():
    random_invalid_press = np.random.uniform(70.1, 100.0)
    max_random_invalid_temp = np.random.uniform(700.1, 1000.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(max_random_invalid_temp, random_invalid_press)


def test_validate_pt_values_negative_range():
    random_nagetive_press = -np.random.uniform(0.1, 100.0)
    random_negative_temp = -np.random.uniform(0.1, 100.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(random_negative_temp, random_nagetive_press)


# -----------------------------------------------------------------------------
