import os
import pathlib
import pickle as pkl

import pytest

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))
DATA_PATH = PATH / "datasets"


@pytest.fixture
def data_path():
    return DATA_PATH.joinpath


@pytest.fixture
def isotherm(data_path):
    path = data_path("isotherm.pkl")
    with open(path, "rb") as f:
        return pkl.load(f)
