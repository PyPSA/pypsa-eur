# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Tests for the split of today's basic chemicals demand into products."""

import pandas as pd
import pytest

from scripts.build_industrial_energy_demand_per_country_today import (
    separate_basic_chemicals,
)

PARAMS = {
    "MWh_NH3_per_tNH3": 5.0,
    "MWh_CH4_per_tNH3_SMR": 10.0,
    "MWh_H2_per_tNH3_electrolysis": 6.0,
    "MWh_elec_per_tNH3_electrolysis": 0.0,
    "MWh_H2_per_tCl": 0.0,
    "MWh_elec_per_tCl": 0.0,
    "MWh_CH4_per_tMeOH": 10.0,
    "MWh_elec_per_tMeOH": 0.0,
}


@pytest.mark.parametrize("endogenous_ammonia", [False, True])
@pytest.mark.parametrize("ammonia, hvc_gas", [(1.0, 20.0), (5.0, 0.0)])
def test_ammonia_gas_removed_from_hvc(endogenous_ammonia, ammonia, hvc_gas):
    carriers = ["gas", "electricity", "hydrogen", "ammonia"]
    index = pd.MultiIndex.from_product([["PL"], carriers])
    demand = pd.DataFrame({"Basic chemicals": [40.0, 5.0, 0.0, 0.0]}, index=index)
    production = pd.DataFrame(
        {"Ammonia": [ammonia], "Methanol": [1.0], "Chlorine": [0.0]}, index=["PL"]
    )

    demand = separate_basic_chemicals(demand, production, PARAMS, endogenous_ammonia)

    assert demand.at[("PL", "gas"), "HVC"] == pytest.approx(hvc_gas)
    assert demand.at[("PL", "gas"), "Ammonia"] == 0.0
