# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Tests for build_electricity_demand.py."""

import numpy as np
import pandas as pd
import pytest

from scripts.build_electricity_demand import reindex_to_snapshots


def hours(year):
    return pd.date_range(f"{year}-01-01", f"{year}-12-31 23:00", freq="h")


@pytest.fixture
def load():
    index = hours(2016).union(hours(2018))
    return pd.DataFrame({"DE": np.arange(len(index), dtype=float)}, index=index)


@pytest.mark.parametrize(
    "snapshot_year, fixed_year",
    [(2018, False), (2018, 2018), (2013, 2018), (2013, 2016), (2016, 2016)],
)
def test_reindex_to_snapshots(load, snapshot_year, fixed_year):
    snapshots = hours(snapshot_year)
    source_year = fixed_year or snapshot_year
    expected = load.loc[snapshots.map(lambda t: t.replace(year=source_year))]

    result = reindex_to_snapshots(load, snapshots, fixed_year)

    assert result.index.equals(snapshots)
    np.testing.assert_array_equal(result.values, expected.values)


def test_reindex_to_snapshots_leap_day_raises(load):
    with pytest.raises(ValueError, match="February 29"):
        reindex_to_snapshots(load, hours(2016), fixed_year=2018)
