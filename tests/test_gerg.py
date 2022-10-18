import numpy as np

import pyforfluids.models as models

import pytest


# -- Components related tests -------------------------------------------------


def test_components():
    model = models.GERG2008()
    valid_components = model.valid_components
    model.validate_components(valid_components)

    with pytest.raises(Exception):
        wrong_component = "other_string"
        model.validate_components(valid_components + [wrong_component])


@pytest.mark.skip("Disabled temporaly")
def test_normalizer():
    model = models.GERG2008()
    composition = {"methane": 0.5, "ethane": 2}

    with pytest.warns(UserWarning):
        model.set_concentration(composition)


# -- Valid ranges tests -------------------------------------------------------


@pytest.mark.skip(reason="no way of currently testing this")
def test_validate_pt_values_no_warnings():
    random = np.random.default_rng(seed=42)
    random_normal_press = random.uniform(0.0, 35.0e6)
    random_normal_temp = random.uniform(90.0, 450.0)
    model = models.GERG2008()
    model.validate_ranges(random_normal_temp, random_normal_press)


def test_validate_pt_values_min_extended_range():
    random = np.random.default_rng(seed=42)
    random_extended_press = random.uniform(35.1e6, 70.0e6)
    min_random_extended_temp = random.uniform(60.0, 89.9)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(min_random_extended_temp, random_extended_press)


def test_validate_pt_values_max_extended_range():
    random = np.random.default_rng(seed=42)
    random_extended_press = random.uniform(35.1e6, 70.0e6)
    max_random_extended_temp = random.uniform(450.1, 700.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(max_random_extended_temp, random_extended_press)


def test_validate_pt_values_min_invalid_range():
    random = np.random.default_rng(seed=42)
    random_invalid_press = random.uniform(70.1e6, 100.0e6)
    min_random_invalid_temp = random.uniform(00.0, 59.9)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(min_random_invalid_temp, random_invalid_press)


def test_validate_pt_values_max_invalid_range():
    random = np.random.default_rng(seed=42)
    random_invalid_press = random.uniform(70.1e6, 100.0e6)
    max_random_invalid_temp = random.uniform(700.1, 1000.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(max_random_invalid_temp, random_invalid_press)


def test_validate_pt_values_negative_range():
    random = np.random.default_rng(seed=42)
    random_nagetive_press = -random.uniform(0.1e6, 100.0e6)
    random_negative_temp = -random.uniform(0.1, 100.0)
    model = models.GERG2008()

    with pytest.warns(UserWarning):
        model.validate_ranges(random_negative_temp, random_nagetive_press)


# -----------------------------------------------------------------------------
