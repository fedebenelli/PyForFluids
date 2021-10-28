import pytest

import numpy as np

import pyforfluids.models as models


@pytest.fixture
def model():
    return models.GERG2008()


# ------------ NORMAL
@pytest.fixture
def random_normal_temp():
    return np.random.uniform(90.0, 450.0)


@pytest.fixture
def random_normal_press():
    return np.random.uniform(0.0, 35.0)


# ------------ EXTENDED
@pytest.fixture
def random_min_extended_temp():
    return np.random.uniform(60.0, 90.0)


@pytest.fixture
def random_max_extended_temp():
    return np.random.uniform(450.0, 700.0)


@pytest.fixture
def random_extended_press():
    return np.random.uniform(35.0, 70.0)


# ----------- INVALID
@pytest.fixture
def random_min_invalid_temp():
    return np.random.uniform(00.0, 60.0)


@pytest.fixture
def random_max_invalid_temp():
    return np.random.uniform(700.0, 1000.0)


@pytest.fixture
def random_invalid_press():
    return np.random.uniform(70.0, 100.0)


# ----------- NEGATIVE
@pytest.fixture
def random_negative_temp():
    return -np.random.uniform(1.0, 100.0)


@pytest.fixture
def random_negative_press():
    return -np.random.uniform(1.0, 100.0)


def test_components(model):
    valid_components = model.valid_components
    model.validate_components(valid_components)

    with pytest.raises(Exception) as e_info:
        wrong_component = "other_string"
        model.validate_components(valid_components + [wrong_component])
    print(e_info)


def test_no_warnings(model, random_normal_press, random_normal_temp):
    model.validate_pt_values(random_normal_press, random_normal_temp)
    assert True


def test_with_warnings_1(
    model, random_extended_press, random_min_extended_temp
):
    with pytest.warns(None) as w_info:
        model.validate_pt_values(
            random_extended_press, random_min_extended_temp
        )
    print(w_info)


def test_with_warnings_2(
    model, random_extended_press, random_max_extended_temp
):
    with pytest.warns(None) as w_info:
        model.validate_pt_values(
            random_extended_press, random_max_extended_temp
        )
    print(w_info)


def test_with_warnings_3(
    model, random_extended_press, random_min_invalid_temp
):
    with pytest.warns(None) as w_info:
        model.validate_pt_values(
            random_extended_press, random_min_invalid_temp
        )
    print(w_info)


def test_with_warnings_4(
    model, random_extended_press, random_max_invalid_temp
):
    with pytest.warns(None) as w_info:
        model.validate_pt_values(
            random_extended_press, random_max_invalid_temp
        )
    print(w_info)


def test_with_warnings_5(model, random_negative_press, random_negative_temp):
    with pytest.warns(None) as w_info:
        model.validate_pt_values(random_negative_press, random_negative_temp)
    print(w_info)
