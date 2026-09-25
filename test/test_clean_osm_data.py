# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the column cleaning functions of scripts/clean_osm_data.py.
"""

import numpy as np
import pandas as pd
import pytest

from scripts.clean_osm_data import (
    _clean_cables,
    _clean_circuits,
    _clean_date,
    _clean_frequency,
    _clean_rating,
    _clean_voltage,
    _clean_wires,
)

# Raw OSM tags arrive as object columns (from JSON) or as pandas' default
# string dtype; both may contain missing values.
DTYPES = [object, "str"]


def _series(values, dtype):
    return pd.Series(values, dtype=dtype)


@pytest.mark.parametrize("dtype", DTYPES)
@pytest.mark.parametrize(
    "func,raw,expected",
    [
        (
            _clean_voltage,
            ["220000", "380 kV", "380000;220000", "medium", np.nan, None],
            ["220000", "380000", "380000;220000", "33000", "", ""],
        ),
        (
            _clean_circuits,
            ["2", "1,5", "1/3", "2;1", np.nan, None],
            ["2", "3", "1", "2;1", "", ""],
        ),
        (
            _clean_cables,
            ["3", "3x2;2", "6;3", "1/3", np.nan, None],
            ["3", "3", "6;3", "1", "", ""],
        ),
        (
            _clean_wires,
            ["single", "Double", "yes", "quad?", np.nan, None],
            ["1", "2", "3", "4", "", ""],
        ),
        (
            _clean_frequency,
            ["50", "50 Hz", "16,7", "16.67", "50;16.7", "0", np.nan, None],
            ["50", "50", "16.7", "16.7", "50;16.7", "0", "", ""],
        ),
    ],
)
def test_clean_string_columns(func, raw, expected, dtype):
    result = func(_series(raw, dtype))
    assert result.tolist() == expected
    assert result.index.equals(pd.RangeIndex(len(raw)))


@pytest.mark.parametrize("dtype", DTYPES)
def test_clean_rating(dtype):
    result = _clean_rating(_series(["1000 MW", "500;500", "700MW"], dtype))
    assert result.tolist() == ["1000", "1000", "700"]


@pytest.mark.parametrize("dtype", DTYPES)
def test_clean_rating_empty_entries(dtype):
    result = _clean_rating(_series(["500;", "unknown", np.nan], dtype))
    assert result.tolist() == ["500", "0", "0"]


@pytest.mark.parametrize("dtype", DTYPES)
def test_clean_date(dtype):
    raw = ["2020-01-01", "circa 2019", "unknown", np.nan, None]
    result = _clean_date(_series(raw, dtype))
    assert pd.api.types.is_datetime64_any_dtype(result)
    assert result.iloc[0] == pd.Timestamp("2020-01-01")
    assert result.iloc[1] == pd.Timestamp("2019-01-01")
    assert result.iloc[2:].isna().all()
