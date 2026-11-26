# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Tests for the processing of technology cost data."""

import pandas as pd
import pytest

from scripts.process_cost_data import scale_solar_cost


@pytest.mark.parametrize(
    ("tech", "dc_ac_ratio", "expected"),
    [
        ("solar-utility", 1.25, 100.0),
        ("solar-utility", 1.5, 116.0),
        ("solar-utility single-axis tracking", 1.0, 84.0),
        ("solar-rooftop", 1.075, 100.0),
    ],
)
def test_scale_solar_cost(tech, dc_ac_ratio, expected):
    index = pd.MultiIndex.from_tuples([(tech, "investment"), (tech, "FOM")])
    costs = pd.DataFrame({"value": [100.0, 2.0]}, index=index)

    costs = scale_solar_cost(costs, tech, dc_ac_ratio, ac_cost_share=0.2)

    assert costs.value.tolist() == pytest.approx([expected, 2.0])
