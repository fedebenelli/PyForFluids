import os
import pathlib
import pickle as pkl

import pandas as pd

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


@pytest.fixture()
def ahmadi323(data_path):
    path = data_path("ahmadi323.3.csv")
    df = pd.read_csv(path)
    df = df.drop("temperature", axis=1)
    return df


@pytest.fixture()
def flash(data_path):
    path = data_path("flash.csv")
    df = pd.read_csv(path)
    df = df.set_index("component")
    return df
