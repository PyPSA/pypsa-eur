# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Tests the functionalities of scripts/build_powerplants.py.
"""

import pandas as pd
import pathlib
import powerplantmatching as pm
import sys

sys.path.append("./scripts")

from test.conftest import get_config_dict
from build_powerplants import add_custom_powerplants, replace_natural_gas_technology, replace_natural_gas_fueltype

path_cwd = pathlib.Path.cwd()


def test_add_custom_powerplants(get_config_dict):
    config_dict = get_config_dict
    config_dict["electricity"]["custom_powerplants"] = True
    custom_powerplants_path = pathlib.Path(path_cwd, "test", "test_data", "custom_powerplants_DE.csv")
    ppl_path = pathlib.Path(path_cwd, "test", "test_data", "powerplants_DE.csv")
    ppl_df = pd.read_csv(ppl_path)
    ppl_final = add_custom_powerplants(ppl_df, custom_powerplants_path, config_dict["electricity"]["custom_powerplants"])
    assert ppl_df.shape == (131, 18)
    assert ppl_final.shape == (137, 18)


