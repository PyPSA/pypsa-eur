# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Tests for the split of today's basic chemicals demand."""

import pandas as pd
import pytest

from scripts.build_industrial_energy_demand_per_country_today import (
    separate_basic_chemicals,
)

PARAMS = {
    "MWh_NH3_per_tNH3": 5.0,
    "MWh_H2_per_tNH3_electrolysis": 6.0,
    "MWh_elec_per_tNH3_electrolysis": 0.25,
    "MWh_CH4_per_tNH3_SMR": 10.0,
    "MWh_elec_per_tNH3_SMR": 0.5,
    "MWh_H2_per_tCl": -1.0,
    "MWh_elec_per_tCl": 3.0,
    "MWh_CH4_per_tMeOH": 10.0,
    "MWh_elec_per_tMeOH": 1.0,
}


@pytest.mark.parametrize("endogenous_ammonia", [True, False])
@pytest.mark.parametrize("gas, hvc_gas", [(40.0, 20.0), (15.0, 0.0)])
def test_ammonia_gas_removed_from_hvc(endogenous_ammonia, gas, hvc_gas):
    carriers = ["gas", "electricity", "hydrogen", "ammonia"]
    index = pd.MultiIndex.from_product([["BE"], carriers])
    demand = pd.DataFrame({"Basic chemicals": [gas, 10.0, 0.0, 0.0]}, index=index)
    production = pd.DataFrame(
        {"Ammonia": 1.0, "Chlorine": 1.0, "Methanol": 1.0}, index=["BE"]
    )

    demand = separate_basic_chemicals(demand, production, PARAMS, endogenous_ammonia)

    hvc = demand.loc["BE", "HVC"]
    assert hvc["gas"] == hvc_gas
    assert hvc["electricity"] == 10.0 - 3.0 - 1.0 - 0.5
    assert hvc["hydrogen"] == 1.0
    ammonia = demand.loc["BE", "Ammonia"]
    assert ammonia["ammonia"] == (5.0 if endogenous_ammonia else 0.0)
    assert ammonia["hydrogen"] == (0.0 if endogenous_ammonia else 6.0)
