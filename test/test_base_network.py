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

from base_network import _get_country, _get_linetypes_config, _get_linetype_by_voltage, _get_oid

path_cwd = pathlib.Path.cwd()


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
        "Al/St 560/50 4-bundle 750.0"
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
        line_type_list.append(
            _get_linetype_by_voltage(v_nom, config["lines"]["types"])
        )
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
