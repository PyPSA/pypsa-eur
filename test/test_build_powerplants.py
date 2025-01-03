# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>>
#
# SPDX-License-Identifier: MIT

"""
Tests the functionalities of scripts/build_powerplants.py.
"""

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.append("./scripts")

from build_powerplants import (
    add_custom_powerplants,
    replace_natural_gas_fueltype,
    replace_natural_gas_technology,
)

path_cwd = pathlib.Path.cwd()


@pytest.mark.parametrize(
    "query_value,expected",
    [(False, (131, 18)), (True, (137, 18))],
)
def test_add_custom_powerplants(config, query_value, expected):
    """
    Verify what returned by add_custom_powerplants.
    """
    config["electricity"]["custom_powerplants"] = query_value
    custom_powerplants_path = pathlib.Path(
        path_cwd, "test", "test_data", "custom_powerplants_DE.csv"
    )
    ppl_path = pathlib.Path(path_cwd, "test", "test_data", "powerplants_DE.csv")
    ppl_df = pd.read_csv(ppl_path)
    ppl_final = add_custom_powerplants(
        ppl_df,
        custom_powerplants_path,
        config["electricity"]["custom_powerplants"],
    )
    assert ppl_df.shape == (131, 18)
    assert ppl_final.shape == expected


def test_replace_natural_gas_technology():
    """
    Verify what returned by replace_natural_gas_technology.
    """
    input_df = pd.DataFrame(
        {
            "Name": [
                "plant_hydro",
                "plant_ng_1",
                "plant_ng_2",
                "plant_ng_3",
                "plant_ng_4",
            ],
            "Fueltype": [
                "Hydro",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
            ],
            "Technology": [
                "Run-Of-River",
                "Steam Turbine",
                "Combustion Engine",
                "Not Found",
                np.nan,
            ],
        }
    )

    reference_df = pd.DataFrame(
        {
            "Name": [
                "plant_hydro",
                "plant_ng_1",
                "plant_ng_2",
                "plant_ng_3",
                "plant_ng_4",
            ],
            "Fueltype": [
                "Hydro",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
                "Natural Gas",
            ],
            "Technology": ["Run-Of-River", "CCGT", "OCGT", "CCGT", "CCGT"],
        }
    )
    modified_df = input_df.assign(Technology=replace_natural_gas_technology)
    comparison_df = modified_df.compare(reference_df)
    assert comparison_df.empty


def test_replace_natural_gas_fueltype():
    """
    Verify what returned by replace_natural_gas_fueltype.
    """
    input_df = pd.DataFrame(
        {
            "Name": [
                "plant_hydro",
                "plant_ng_1",
                "plant_ng_2",
            ],
            "Fueltype": [
                "Hydro",
                "Gas",
                "Natural",
            ],
            "Technology": [
                "Run-Of-River",
                "CCGT",
                "OCGT",
            ],
        }
    )

    reference_df = pd.DataFrame(
        {
            "Name": [
                "plant_hydro",
                "plant_ng_1",
                "plant_ng_2",
            ],
            "Fueltype": [
                "Hydro",
                "Natural Gas",
                "Natural Gas",
            ],
            "Technology": [
                "Run-Of-River",
                "CCGT",
                "OCGT",
            ],
        }
    )
    modified_df = input_df.assign(Fueltype=replace_natural_gas_fueltype)
    comparison_df = modified_df.compare(reference_df)
    assert comparison_df.empty
