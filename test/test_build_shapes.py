# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the functionalities of scripts/build_shapes.py.
"""

import pathlib
import sys

import geopandas as gpd
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
                301184483954.051,
                "POINT (1170973.133650994 4959324.697171187)",
                7715732.183141501,
            ),
        ),
        (
            10.0,
            (
                301184417650.59454,
                "POINT (1170973.0551039157 4959324.750552279)",
                7715728.487943892,
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
    gdf_country_simplified["centroid"] = gdf_country_simplified.centroid
    gdf_country_simplified["perimeter"] = gdf_country_simplified.length
    output_tuple = (
        gdf_country_simplified["area"][0],
        str(gdf_country_simplified["centroid"][0]),
        gdf_country_simplified["perimeter"][0],
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
                89.48805178772852,
                "POINT (12.612315614274285 40.81446088016855)",
                61.66613129515858,
            ),
        ),
        (
            ["DE"],
            (
                53.71978281859542,
                "POINT (10.08343387158037 51.60752846577776)",
                56.456408985211766,
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
    europe_shape_gdf["centroid"] = europe_shape_gdf.centroid
    europe_shape_gdf["perimeter"] = europe_shape_gdf.length
    output_tuple = (
        europe_shape_gdf["area"][0],
        str(europe_shape_gdf["centroid"][0]),
        europe_shape_gdf["perimeter"][0],
    )
    print(output_tuple)
    assert len(output_tuple) == len(expected_tuple)
    assert all([x == y for x, y in zip(output_tuple, expected_tuple)])
