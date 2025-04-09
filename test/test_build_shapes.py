# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the functionalities of scripts/build_shapes.py.
"""

import pathlib
import sys

import geopandas as gpd
import numpy as np
import pytest

sys.path.append("./scripts")

from build_shapes import _simplify_polys, eez

path_cwd = pathlib.Path.cwd()


@pytest.mark.parametrize(
    "tolerance,expected_tuple",
    [
        (
            None,
            (
                301184483954.05,
                7715732.18,
            ),
        ),
        (
            10.0,
            (
                301184417650.59,
                7715728.49,
            ),
        ),
    ],
)
def test_simplify_polys(tolerance, expected_tuple, italy_shape):
    """
    Verify what is returned by _simplify_polys.

    Note:
        - tolerance = None, no simplification takes place
    """
    gdf_country = gpd.read_file(italy_shape).to_crs(6933)
    simplified_polys = _simplify_polys(gdf_country.geometry, tolerance=tolerance)
    gdf_country_simplified = gpd.GeoDataFrame(geometry=simplified_polys)
    gdf_country_simplified["area"] = gdf_country_simplified.area
    gdf_country_simplified["perimeter"] = gdf_country_simplified.length
    output_tuple = (
        np.round(gdf_country_simplified["area"][0], 2),
        np.round(gdf_country_simplified["perimeter"][0], 2),
    )
    assert len(output_tuple) == len(expected_tuple)
    assert all([x == y for x, y in zip(output_tuple, expected_tuple)])


@pytest.mark.parametrize(
    "country_list",
    [["DE"], ["IT"]],
)
def test_eez(config, country_list, download_eez):
    """
    Verify what is returned by eez.
    """
    eez_path = download_eez
    offshore_shapes_gdf = eez(eez_path, country_list)
    assert offshore_shapes_gdf.shape == (1, 1)
    assert offshore_shapes_gdf.index == country_list[0]
