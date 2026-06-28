# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Regression tests for build_electricity_demand.py.

Tests for the fixed_year reindexing logic that was silently producing
NaN load data when fixed_year differed from the snapshot year.
See https://github.com/PyPSA/pypsa-eur/issues/2187
"""

import numpy as np
import pandas as pd
import pytest

from scripts.build_electricity_demand import reindex_load_fixed_year


def _make_load(year, countries=("DE", "FR")):
    """Create synthetic hourly load data for a given year."""
    idx = pd.date_range(f"{year}-01-01", f"{year}-12-31 23:00", freq="h")
    rng = np.random.default_rng(42)
    data = {c: rng.random(len(idx)) * 1000 + 500 for c in countries}
    return pd.DataFrame(data, index=idx)


class TestFixedYearReindex:
    """Regression tests for #2187: fixed_year must not produce NaN load."""

    def test_same_year(self):
        """fixed_year == snapshot year: trivial case, must work."""
        load = _make_load(2018)
        snapshots = pd.date_range("2018-01-01", "2018-12-31 23:00", freq="h")
        result = reindex_load_fixed_year(load, snapshots, fixed_year=2018)
        assert not result.isna().any().any()
        assert len(result) == len(snapshots)

    def test_different_year_no_nans(self):
        """Core regression test for #2187: fixed_year != snapshot year."""
        load = _make_load(2018)
        snapshots = pd.date_range("2013-01-01", "2013-12-31 23:00", freq="h")
        result = reindex_load_fixed_year(load, snapshots, fixed_year=2018)
        assert not result.isna().any().any(), "fixed_year produced NaN load data"
        assert len(result) == len(snapshots)

    def test_leap_fixed_year_to_non_leap_snapshots(self):
        """Leap fixed_year (2024) mapped to non-leap snapshots (2013)."""
        load = _make_load(2024)
        snapshots = pd.date_range("2013-01-01", "2013-12-31 23:00", freq="h")
        result = reindex_load_fixed_year(load, snapshots, fixed_year=2024)
        assert not result.isna().any().any()
        assert len(result) == 8760  # non-leap year hours

    def test_non_leap_fixed_year_with_leap_snapshots_raises(self):
        """Non-leap fixed_year (2018) with leap snapshots containing Feb 29."""
        load = _make_load(2018)
        # Leap year snapshots WITH Feb 29
        snapshots = pd.date_range("2024-01-01", "2024-12-31 23:00", freq="h")
        with pytest.raises(ValueError, match="Feb 29"):
            reindex_load_fixed_year(load, snapshots, fixed_year=2018)

    def test_no_fixed_year(self):
        """No fixed_year: standard reindex, must not regress."""
        load = _make_load(2018)
        snapshots = pd.date_range("2018-01-01", "2018-12-31 23:00", freq="h")
        result = reindex_load_fixed_year(load, snapshots, fixed_year=False)
        assert not result.isna().any().any()

    def test_values_are_preserved(self):
        """Verify actual load values are mapped correctly, not just non-NaN."""
        load = _make_load(2018)
        snapshots = pd.date_range("2013-01-01", "2013-12-31 23:00", freq="h")
        result = reindex_load_fixed_year(load, snapshots, fixed_year=2018)
        # First hour of 2013 should have the value from first hour of 2018
        assert result.iloc[0]["DE"] == load.iloc[0]["DE"]
        # Last hour likewise
        assert result.iloc[-1]["FR"] == load.iloc[-1]["FR"]
