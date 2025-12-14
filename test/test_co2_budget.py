# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import pytest

from scripts.co2_budget import bound_value_for_horizon, co2_budget_for_horizon


def test_bound_value_for_horizon_null():
    assert bound_value_for_horizon(None, 2030) is None


def test_bound_value_for_horizon_scalar():
    assert bound_value_for_horizon(0.1, 2030) == 0.1


def test_bound_value_for_horizon_dict_year_int():
    bound = {2030: 2.0}
    assert bound_value_for_horizon(bound, 2030) == 2.0
    assert bound_value_for_horizon(bound, 2040) is None


def test_bound_value_for_horizon_dict_year_str():
    bound = {"2030": 2.0}
    assert bound_value_for_horizon(bound, 2030) == 2.0


def test_co2_budget_for_horizon_absolute_dict_and_null_lower():
    co2_budget = {
        "values": "absolute",
        "upper": {2050: 0.1},
        "lower": None,
    }
    upper, lower = co2_budget_for_horizon(co2_budget, current_horizon=2050)
    assert upper == 0.1
    assert lower is None


def test_co2_budget_for_horizon_fraction_requires_baseline():
    co2_budget = {
        "values": "fraction",
        "upper": 0.5,
        "lower": None,
    }
    with pytest.raises(ValueError, match="1990 baseline"):
        co2_budget_for_horizon(co2_budget, current_horizon=2030)


def test_co2_budget_for_horizon_fraction_applies_baseline():
    co2_budget = {
        "values": "fraction",
        "upper": 0.5,
        "lower": 0.1,
    }
    upper, lower = co2_budget_for_horizon(
        co2_budget, current_horizon=2030, baseline_1990=10.0
    )
    assert upper == 5.0
    assert lower == 1.0


def test_co2_budget_for_horizon_rejects_lower_without_upper():
    co2_budget = {
        "values": "absolute",
        "upper": None,
        "lower": 0.1,
    }
    with pytest.raises(ValueError, match="requires an upper constraint"):
        co2_budget_for_horizon(co2_budget, current_horizon=2030)


def test_co2_budget_for_horizon_rejects_lower_greater_equal_upper():
    co2_budget = {
        "values": "absolute",
        "upper": 0.1,
        "lower": 0.1,
    }
    with pytest.raises(ValueError, match="Lower bound"):
        co2_budget_for_horizon(co2_budget, current_horizon=2030)
