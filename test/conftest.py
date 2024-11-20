# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8

import pathlib
import shutil

import pypsa
import pytest
import yaml


@pytest.fixture(scope="function")
def get_power_network_scigrid_de():
    return pypsa.examples.scigrid_de(from_master=True)


@pytest.fixture(scope="function")
def get_power_network_ac_dc_meshed():
    return pypsa.examples.ac_dc_meshed(from_master=True)


@pytest.fixture(scope="function")
def get_config_dict():
    path_config = pathlib.Path(pathlib.Path.cwd(), "config", "config.default.yaml")
    with open(path_config, "r") as file:
        config_dict = yaml.safe_load(file)
    return config_dict
