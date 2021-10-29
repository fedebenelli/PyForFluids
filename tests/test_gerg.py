import numpy as np

import pyforfluids.models as models

import pytest


def test_components():
    model = models.GERG2008()
    valid_components = model.valid_components
    model.validate_components(valid_components)

    with pytest.raises(Exception) as e_info:
        wrong_component = "other_string"
        model.validate_components(valid_components + [wrong_component])
    print(e_info)


def test_validate_pt_values_no_warnings():
    random_normal_press = np.random.uniform(0.0, 35.0)
    random_normal_temp = np.random.uniform(90.0, 450.0)
    model = models.GERG2008()
    model.validate_pt_values(random_normal_press, random_normal_temp)
    assert True


def test_validate_pt_values_min_extended_range():
    random_extended_press = np.random.uniform(35.1, 70.0)
    min_random_extended_temp = np.random.uniform(60.0, 89.9)
    model = models.GERG2008()
    with pytest.warns(UserWarning) as w_info:
        model.validate_pt_values(
            random_extended_press, min_random_extended_temp
        )
    print(w_info)


def test_validate_pt_values_max_extended_range():
    random_extended_press = np.random.uniform(35.1, 70.0)
    max_random_extended_temp = np.random.uniform(450.1, 700.0)
    model = models.GERG2008()
    with pytest.warns(UserWarning) as w_info:
        model.validate_pt_values(
            random_extended_press, max_random_extended_temp
        )
    print(w_info)


def test_validate_pt_values_min_invalid_range():
    random_invalid_press = np.random.uniform(70.1, 100.0)
    min_random_invalid_temp = np.random.uniform(00.0, 59.9)
    model = models.GERG2008()
    with pytest.warns(UserWarning) as w_info:
        model.validate_pt_values(random_invalid_press, min_random_invalid_temp)
    print(w_info)


def test_validate_pt_values_max_invalid_range():
    random_invalid_press = np.random.uniform(70.1, 100.0)
    max_random_invalid_temp = np.random.uniform(700.1, 1000.0)
    model = models.GERG2008()
    with pytest.warns(UserWarning) as w_info:
        model.validate_pt_values(random_invalid_press, max_random_invalid_temp)
    print(w_info)


def test_validate_pt_values_negative_range():
    random_nagetive_press = -np.random.uniform(0.1, 100.0)
    random_negative_temp = -np.random.uniform(0.1, 100.0)
    model = models.GERG2008()
    with pytest.warns(UserWarning) as w_info:
        model.validate_pt_values(random_nagetive_press, random_negative_temp)
    print(w_info)
