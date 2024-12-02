# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8

import pandas as pd
import pathlib

import pypsa
import pytest
import yaml


@pytest.fixture(scope="function")
def scigrid_network():
    return pypsa.examples.scigrid_de(from_master=True)


@pytest.fixture(scope="function")
def ac_dc_network():
    return pypsa.examples.ac_dc_meshed(from_master=True)


@pytest.fixture(scope="function")
def config():
    path_config = pathlib.Path(pathlib.Path.cwd(), "config", "config.default.yaml")
    with open(path_config, "r") as file:
        config_dict = yaml.safe_load(file)
    return config_dict


@pytest.fixture(scope="function")
def converters_dataframe():
    return pd.DataFrame(
        {
            "index": 0,
            "converter_id": "convert_20_41",
            "bus0": "41",
            "bus1": "42",
            "underground": False,
            "under_construction": False,
            "country": "US",
            "geometry": "LINESTRING(-122.3787 37.6821, -122.3777 37.6831)",
        },
        index=[0],
    ).set_index("converter_id")


@pytest.fixture(scope="function")
def buses_dataframe():
    return pd.DataFrame(
        {
            "bus_id": [5231, 5232],
            "voltage": [380.0, 400.0],
            "dc": ["f", "f"],
            "symbol": ["Substation", "Substation"],
            "under_construction": ["f", "f"],
            "x": [6.8884, 6.8894],
            "y": [45.6783, 45.6793],
            "country": ["IT", "IT"],
            "geometry": ["POINT (6.8884 45.6783)", "POINT (6.8894 45.6793)"],
        },
        index=[5090, 5091],
    )
