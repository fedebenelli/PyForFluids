import os
import pathlib
import pickle as pkl
from collections import OrderedDict

import numpy as np

import pandas as pd

from pyforfluids.models import CubicEOS

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
def peng_robinson_pv(data_path):
    path = data_path("peng_robinson_mix_pv.csv")
    df = pd.read_csv(path)
    return df


@pytest.fixture()
def peng_robinson_mix():
    components = [
        "CO2",
        "C1-N2",
        "C2-C3",
        "C4",
        "C5",
        "C6",
        "C7+",
        "Asph",
    ]

    z = np.array(
        [0.0246, 0.3694, 0.0752, 0.0193, 0.0157, 0.0162, 0.4633, 0.0163],
        dtype="d",
    )

    composition = OrderedDict()
    for name, conc in zip(components, z):
        composition[name] = conc

    tc = np.array(
        [
            304.039,
            189.428,
            339.872,
            419.817,
            465.094,
            507.317,
            860.372,
            1424.817,
        ],
        dtype="d",
    )
    pc = np.array(
        [73.790, 45.830, 45.410, 37.540, 33.800, 32.900, 12.460, 12.290],
        dtype="d",
    )
    w = np.array(
        [
            0.22500000000000001,
            0.0085000000000000006,
            0.12709999999999999,
            0.18779999999999999,
            0.2397,
            0.27500000000000002,
            1.022,
            1.4410000000000001,
        ],
        dtype="d",
    )
    kij = np.array(
        [
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 5.23e-002, 0.135],
            [0, 0, 0, 0, 0, 0, 0, 0.135],
            [0, 0, 0, 0, 0, 0, 0, 0.135],
            [0, 0, 0, 0, 0, 0, 0, 0.135],
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, 5.3e-002, 0, 0, 0, 0, 0, 0],
            [0, 0.135, 0.135, 0.135, 0.135, 0, 0, 0],
        ],
        dtype="d",
        order="F",
    )
    pr = CubicEOS(
        model="PR",
        mix_rule="ClassicVdW",
        names=components,
        critical_temperature=tc,
        critical_pressure=pc,
        acentric_factor=w,
        kij_matrix=kij,
        lij_matrix=0 * kij,
    )
    return composition, pr


@pytest.fixture()
def flash(data_path):
    path = data_path("flash.csv")
    df = pd.read_csv(path)
    df = df.set_index("component")
    return df
