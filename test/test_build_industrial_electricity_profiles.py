# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Tests for build_industrial_electricity_profiles.py."""

import numpy as np
import pandas as pd
import pytest

from scripts.build_industrial_electricity_profiles import (
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


@pytest.mark.parametrize("year", [2013, 2020])
def test_shift_profiles_aligns_weekdays(ffe_profiles, year):
    shifted = shift_profiles(ffe_profiles, year, "GB").iloc[:, 0]

    assert shifted.index.equals(hours(year))
    assert np.isclose(shifted.mean(), 1)

    weekday_levels = shifted.groupby(shifted.index.dayofweek).median()
    assert np.allclose(weekday_levels / weekday_levels[0], np.arange(1, 8))


def test_build_electricity_profiles(ffe_profiles):
    nodes = ["DE0 0", "FI0 0", "GB0 0"]
    data = [[3.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
    demand = pd.DataFrame(data, index=["Cement", "Ammonia"], columns=nodes)
    week = pd.date_range("2013-01-14", "2013-01-20 23:00", freq="h")

    full_year = build_electricity_profiles(demand, ffe_profiles, hours(2013))
    profiles = build_electricity_profiles(demand, ffe_profiles, week)

    # annual mean load of a flat profile, also without demand
    assert np.allclose(full_year.mean() * 8760 / 1e6, demand.sum())
    # partial years keep the level of the full year
    pd.testing.assert_frame_equal(profiles, full_year.loc[week])
    # noon peak in local time
    assert list(profiles.iloc[:, :2].idxmax().dt.hour) == [11, 10]
