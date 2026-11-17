# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Tests for the split of industry fuel demand into temperature bands."""

import pandas as pd
import pytest

from scripts.build_industry_sector_ratios_endogenous import (
    FUELS,
    build_industry_sector_ratios_endogenous,
)

BANDS = ["heat", "heat100-200", "heat200-500", "heat500+"]


@pytest.fixture
def ratios():
    sectors = ["Cement", "HVC", "Electric arc"]
    columns = pd.MultiIndex.from_product([["BE"], sectors])
    demand = pd.DataFrame(2.0, index=FUELS + ["electricity"], columns=columns)
    production = pd.DataFrame(1e3, index=["BE"], columns=sectors + ["Ammonia"])
    shares = pd.DataFrame(
        [[0.1, 0.2, 0.3, 0.4], [0.5, 0.5, 0.0, 0.0]],
        index=["Cement", "HVC"],
        columns=["heat<100", "heat100-200", "heat200-500", "heat500+"],
    )
    return build_industry_sector_ratios_endogenous(demand, production, shares, 1.0)[
        "BE"
    ]


@pytest.mark.parametrize("sector, heat", [("Cement", 14.0), ("HVC", 11.0)])
def test_fuels_become_heat(ratios, sector, heat):
    assert list(ratios.columns) == ["Cement", "HVC"]
    assert ratios.loc[BANDS, sector].sum() == pytest.approx(heat)
    assert ratios.at["elec", sector] == 2.0
    assert ratios.loc[["coke", "biomass", "methane"], sector].sum() == 0.0
    assert pd.isna(ratios.at["naphtha", sector]) == (sector == "HVC")
