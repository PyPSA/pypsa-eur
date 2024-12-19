# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Tests the functionalities of scripts/build_shapes.py.
"""

import pathlib
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest

sys.path.append("./scripts")

from build_shapes import (
    _simplify_polys
)

path_cwd = pathlib.Path.cwd()


@pytest.mark.parametrize(
    "tolerance,expected_tuple",
    [
        (None, (837421026967.6136, "POINT (1222256.2540812986 4769376.881403567)", 5958197.627333247)),
        (10.0, (837420795181.8837, "POINT (1222256.3488327404 4769376.890839249)", 5958195.650234007))],
)
def test_simplify_polys(tolerance, expected_tuple):
    """
    Verify what is returned by _simplify_polys.

    Note:
        - tolerance = None, no simplification takes place
    """
    italy_shape = pathlib.Path(path_cwd, "test", "test_data", "italy_shape.geojson")
    gdf_country = gpd.read_file(italy_shape).to_crs(6933)
    simplified_polys = _simplify_polys(gdf_country.geometry, tolerance=tolerance)
    gdf_country_simplified = gpd.GeoDataFrame(geometry=simplified_polys)
    gdf_country_simplified["area"] = gdf_country_simplified.area
    gdf_country_simplified["centroid"] = gdf_country_simplified.centroid
    gdf_country_simplified["perimeter"] = gdf_country_simplified.length
    output_tuple = (gdf_country_simplified["area"][0], str(gdf_country_simplified["centroid"][0]), gdf_country_simplified["perimeter"][0])
    print(output_tuple)
    assert len(output_tuple) == len(expected_tuple)
    assert all([x == y for x, y in zip(output_tuple, expected_tuple)])
