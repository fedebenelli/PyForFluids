import pytest

import pyforfluids.models as models


@pytest.fixture
def model():
    return models.GERG2008()


def test_components(model):
    valid_components = model.valid_components
    model.validate_components(valid_components)

    with pytest.raises(Exception) as e_info:
        wrong_component = "other_string"
        model.validate_components(valid_components + [wrong_component])
    print(e_info)
