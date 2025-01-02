# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8

import pathlib

import pandas as pd
import pypsa
import pytest
import yaml


@pytest.fixture(scope="function")
def scigrid_network():
    return pypsa.examples.scigrid_de(from_master=True)


@pytest.fixture(scope="function")
def ac_dc_network():
    return pypsa.examples.ac_dc_meshed(from_master=True)


@pytest.fixture(scope="session")
def config():
    path_config = pathlib.Path(pathlib.Path.cwd(), "config", "config.default.yaml")
    with open(path_config, "r") as file:
        config_dict = yaml.safe_load(file)
    return config_dict


@pytest.fixture(scope="function")
def buses_dataframe():
    return pd.DataFrame(
        {
            "bus_id": [5231, 5232],
            "voltage": [380.0, 380.0],
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


@pytest.fixture(scope="function")
def converters_dataframe():
    return pd.DataFrame(
        {
            "converter_id": "convert_5231_5232",
            "bus0": 5231,
            "bus1": 5232,
            "voltage": 380.0,
            "geometry": "'LINESTRING(6.8884 45.6783 ",
            "": "6.8894 45.6793)'",
        },
        index=[0],
    )


@pytest.fixture(scope="function")
def lines_dataframe():
    return pd.DataFrame(
        {
            "line_id": "line_5231_5232",
            "bus0": 5231,
            "bus1": 5232,
            "voltage": 380.0,
            "circuits": 1.0,
            "length": 1000.0,
            "underground": "t",
            "under_construction": "f",
            "geometry": "'LINESTRING(6.8884 45.6783 ",
            "": "6.8894 45.6793)'",
        },
        index=[0],
    )


@pytest.fixture(scope="function")
def links_dataframe():
    return pd.DataFrame(
        {
            "link_id": "link_5231_5232",
            "bus0": 5231,
            "bus1": 5232,
            "voltage": 380.0,
            "p_nom": 600.0,
            "length": 1000.0,
            "underground": "t",
            "under_construction": "f",
            "geometry": "'LINESTRING(6.8884 45.6783 ",
            "": "6.8894 45.6793)'",
        },
        index=[0],
    )


@pytest.fixture(scope="function")
def transformers_dataframe():
    return pd.DataFrame(
        {
            "transformer_id": "transf_5231_5232",
            "bus0": 5231,
            "bus1": 5232,
            "voltage": 380.0,
            "geometry": "'LINESTRING(6.8884 45.6783 ",
            "": "6.8894 45.6793)'",
        },
        index=[0],
    )
