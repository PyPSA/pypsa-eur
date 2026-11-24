# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Tests for build_industrial_energy_demand_per_node.py."""

import numpy as np
import pandas as pd
import pytest

from scripts.build_industrial_energy_demand_per_node import load_sector_ratios

SECTORS = ["Cement", "HVC", "Electric arc"]


@pytest.fixture
def files(tmp_path):
    columns = pd.MultiIndex.from_product([["BE"], SECTORS])
    horizon = pd.DataFrame(
        [[3.0] * 3, [2.0] * 3, [1.0] * 3],
        index=["elec", "methane", "naphtha"],
        columns=columns,
    )
    endogenous = pd.DataFrame(
        [[1.0, 1.0], [0.0, 0.0], [0.0, np.nan], [4.0, 3.0]],
        index=["elec", "methane", "naphtha", "heat100-200"],
        columns=columns[:2],
    )
    horizon.to_csv(tmp_path / "horizon.csv")
    endogenous.to_csv(tmp_path / "endogenous.csv")
    return tmp_path / "horizon.csv", tmp_path / "endogenous.csv"


@pytest.mark.parametrize("industry_heat", [True, False])
def test_load_sector_ratios(files, industry_heat):
    fn, fn_endogenous = files
    ratios = load_sector_ratios(fn, fn_endogenous if industry_heat else [], ["cement"])[
        "BE"
    ]

    if industry_heat:
        expected = [[0.0, 1.0, 3.0], [0.0, 0.0, 2.0], [0.0, 1.0, 1.0], [0.0, 3.0, 0.0]]
        index = ["elec", "methane", "naphtha", "heat100-200"]
    else:
        expected = [[0.0, 3.0, 3.0], [0.0, 2.0, 2.0], [0.0, 1.0, 1.0]]
        index = ["elec", "methane", "naphtha"]
    expected = pd.DataFrame(expected, index=index, columns=SECTORS)

    pd.testing.assert_frame_equal(ratios, expected, check_names=False)
