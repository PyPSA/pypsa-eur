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
    _get_oid,
    _get_country
)

path_cwd = pathlib.Path.cwd()


@pytest.mark.parametrize(
    "column_name, expected",
    [("tags", pd.Series(["12", "12345", "9876654"])), ("other", pd.Series([np.nan, np.nan, np.nan]))],
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


@pytest.mark.parametrize(
    "column_name, expected",
    [("tags", pd.Series(["NG", "CH", "AU"])), ("other", pd.Series([np.nan, np.nan, np.nan]))],
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
