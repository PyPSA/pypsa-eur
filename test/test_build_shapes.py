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

from build_shapes import _simplify_polys, countries, country_cover, eez

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
    [["MK"], ["IT"]],
)
def test_countries(config, download_natural_earth, country_list):
    """
    Verify what is returned by countries.
    """
    natural_earth = download_natural_earth
    country_shapes_df = countries(natural_earth, country_list)
    assert country_shapes_df.shape == (1,)
    assert country_shapes_df.index.unique().tolist() == country_list


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


@pytest.mark.parametrize(
    "country_list,expected_tuple",
    [
        (
            ["IT"],
            (
                89.49,
                61.67,
            ),
        ),
        (
            ["DE"],
            (
                53.72,
                56.46,
            ),
        ),
    ],
)
def test_country_cover(
    country_list, download_natural_earth, download_eez, expected_tuple
):
    """
    Verify what is returned by country_cover.
    """
    natural_earth = download_natural_earth
    eez_path = download_eez
    country_shapes_gdf = countries(natural_earth, country_list).reset_index()
    offshore_shapes_gdf = eez(eez_path, country_list).reset_index()
    europe_shape_gdf = gpd.GeoDataFrame(
        geometry=[country_cover(country_shapes_gdf, offshore_shapes_gdf.geometry)],
        crs=6933,
    )
    europe_shape_gdf["area"] = europe_shape_gdf.area
    europe_shape_gdf["perimeter"] = europe_shape_gdf.length
    output_tuple = (
        np.round(europe_shape_gdf["area"][0], 2),
        np.round(europe_shape_gdf["perimeter"][0], 2),
    )
    assert len(output_tuple) == len(expected_tuple)
    assert all([x == y for x, y in zip(output_tuple, expected_tuple)])
