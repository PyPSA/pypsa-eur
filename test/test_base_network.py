# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Tests the functionalities of scripts/base_network.py.
"""

import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.append("./scripts")

from base_network import (
    _get_country,
    _get_linetype_by_voltage,
    _get_linetypes_config,
    _get_oid,
    _load_buses,
    _load_converters_from_osm,
)

path_cwd = pathlib.Path.cwd()


df_converters_reference = pd.DataFrame(
    {
        "converter_id": "convert_20_41",
        "index": 0,
        "bus0": "41",
        "bus1": "42",
        "underground": False,
        "under_construction": False,
        "country": "US",
        "geometry": "LINESTRING(-122.3787 37.6821, -122.3777 37.6831)",
        "carrier": "B2B",
        "dc": True,
    },
    index=[0],
).set_index("converter_id")


@pytest.mark.parametrize(
    "column_name, expected",
    [
        ("tags", pd.Series(["NG", "CH", "AU"])),
        ("other", pd.Series([np.nan, np.nan, np.nan])),
    ],
)
def test_get_country(column_name, expected):
    """
    Verify what returned by _get_country()
    """
    data_list = [['"country"=>"NG"'], ['"country"=>"CH"'], ['"country"=>"AU"']]
    df_exercise = pd.DataFrame(data_list, columns=[column_name])
    output_series = _get_country(df_exercise)
    comparison_series = output_series.compare(expected)
    assert comparison_series.size == 0


def test_get_linetypes_config(config):
    """
    Verify what returned by _get_linetypes_config.
    """
    output_dict = _get_linetypes_config(
        config["lines"]["types"], config["electricity"]["voltages"]
    )
    assert output_dict == config["lines"]["types"]


def test_get_linetype_by_voltage(config):
    """
    Verify what returned by _get_linetype_by_voltage.
    """

    reference_list = [
        "Al/St 240/40 2-bundle 220.0",
        "Al/St 240/40 3-bundle 300.0",
        "Al/St 240/40 3-bundle 300.0",
        "Al/St 240/40 4-bundle 380.0",
        "Al/St 240/40 4-bundle 380.0",
        "Al/St 240/40 4-bundle 380.0",
        "Al/St 560/50 4-bundle 750.0",
    ]

    v_nom_list = [
        220.0,
        300.0,
        330.0,
        380.0,
        400.0,
        500.0,
        750.0,
    ]

    line_type_list = []

    for v_nom in v_nom_list:
        line_type_list.append(_get_linetype_by_voltage(v_nom, config["lines"]["types"]))
    assert len(line_type_list) == len(reference_list)
    assert all([x == y for x, y in zip(line_type_list, reference_list)])


@pytest.mark.parametrize(
    "column_name, expected",
    [
        ("tags", pd.Series(["12", "12345", "9876654"])),
        ("other", pd.Series([np.nan, np.nan, np.nan])),
    ],
)
def test_get_oid(column_name, expected):
    """
    Verify what returned by _get_oid()
    """
    data_list = [['"oid"=>"12"'], ['"oid"=>"12345"'], ['"oid"=>"9876654"']]
    df_exercise = pd.DataFrame(data_list, columns=[column_name])
    output_series = _get_oid(df_exercise)
    print(output_series)
    comparison_series = output_series.compare(expected)
    assert comparison_series.size == 0


# def test_load_converters_from_osm(tmpdir, buses_dataframe, config, converters_dataframe):
#     """
#     Verify what returned by _load_converters_from_osm.
#     """
#     buses_path = pathlib.Path(tmpdir, "buses.csv")
#     buses_dataframe.to_csv(buses_path, index=False)
#     countries = config["countries"]
#     europe_shape = pathlib.Path(path_cwd, "resources", "europe_shape.geojson")
#     df_buses = _load_buses(buses_path, europe_shape, countries, config).reset_index()
#     converters_path = pathlib.Path(tmpdir, "converters_exercise.csv")
#     converters_dataframe.to_csv(converters_path, index=False)
#     df_converters_output = _load_converters_from_osm(df_buses, converters_path)
#     df_converters_comparison = df_converters_output.compare(df_converters_reference)
#     pathlib.Path.unlink(df_converters_comparison)
#     assert df_converters_comparison.empty


def test_load_buses(tmpdir, config, buses_dataframe):
    """
    Verify what returned by _load_buses.
    """
    df_buses_reference = pd.DataFrame(
        {
            "bus_id": ["5231", "5232"],
            "v_nom": [380.0, 400.0],
            "symbol": ["Substation", "Substation"],
            "under_construction": [False, False],
            "x": [6.8884, 6.8894],
            "y": [45.6783, 45.6793],
            "country": ["IT", "IT"],
            "geometry": ["POINT (6.8884 45.6783)", "POINT (6.8894 45.6793)"],
            "carrier": ["AC", "AC"],
        },
    )
    buses_path = pathlib.Path(tmpdir, "buses.csv")
    buses_dataframe.to_csv(buses_path, index=False)
    countries = config["countries"]
    europe_shape = pathlib.Path(path_cwd, "resources", "europe_shape.geojson")
    df_buses_output = _load_buses(
        buses_path, europe_shape, countries, config
    ).reset_index()
    pathlib.Path.unlink(buses_path)
    df_comparison = df_buses_output.compare(df_buses_reference)
    assert df_comparison.empty
