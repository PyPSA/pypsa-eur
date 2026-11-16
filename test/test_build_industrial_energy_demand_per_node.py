# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Tests for build_industrial_energy_demand_per_node.py."""

import numpy as np
import pandas as pd
import pytest

from scripts.build_industrial_energy_demand_per_node import (
    FFE_PROFILES,
    FFE_YEAR,
    build_electricity_profiles,
    shift_profiles,
)


def hours(year):
    return pd.date_range(f"{year}-01-01", f"{year}-12-31 23:00", freq="h")


@pytest.fixture
def ffe_profiles():
    """Profiles with a weekday level and a peak at noon local time."""
    index = hours(FFE_YEAR)
    values = 1.0 + index.dayofweek + (index.hour == 12)
    return pd.DataFrame(dict.fromkeys(FFE_PROFILES.values(), values), index=index)


@pytest.fixture
def demand():
    nodes = ["DE0 0", "FI0 0", "GB0 0", "GB0 1"]
    data = [[3.0, 1.0, 2.0, 0.0], [1.0, 0.0, 4.0, 0.0]]
    return pd.DataFrame(data, index=["Cement", "Ammonia"], columns=nodes)


@pytest.mark.parametrize("year", [2013, 2019, 2020, 2024])
@pytest.mark.parametrize("country", ["DE", "GB"])
def test_shift_profiles(ffe_profiles, year, country):
    shifted = shift_profiles(ffe_profiles, year, country).iloc[:, 0]

    assert shifted.index.equals(hours(year))
    assert np.isclose(shifted.mean(), 1)

    weekday_levels = shifted.groupby(shifted.index.dayofweek).median()
    assert np.allclose(weekday_levels / weekday_levels[0], np.arange(1, 8))


@pytest.mark.parametrize("year", [2013, 2020])
def test_build_electricity_profiles_keeps_mean_load(ffe_profiles, demand, year):
    profiles = build_electricity_profiles(demand, ffe_profiles, hours(year))

    assert profiles.index.equals(hours(year))
    assert np.allclose(profiles.mean() * 8760 / 1e6, demand.sum(), rtol=1e-3)


def test_build_electricity_profiles_partial_year(ffe_profiles, demand):
    full_year = build_electricity_profiles(demand, ffe_profiles, hours(2013))
    week = pd.date_range("2013-03-01", "2013-03-07 23:00", freq="h")

    profiles = build_electricity_profiles(demand, ffe_profiles, week)

    pd.testing.assert_frame_equal(profiles, full_year.loc[week])


def test_build_electricity_profiles_local_time(ffe_profiles, demand):
    day = pd.date_range("2013-01-15", "2013-01-15 23:00", freq="h")

    profiles = build_electricity_profiles(demand, ffe_profiles, day)

    assert list(profiles.iloc[:, :3].idxmax().dt.hour) == [11, 10, 12]


def test_build_electricity_profiles_unknown_subsector(ffe_profiles, demand):
    with pytest.raises(KeyError):
        build_electricity_profiles(
            demand.rename(index={"Cement": "Unknown"}), ffe_profiles, hours(2013)
        )
