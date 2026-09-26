# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Tests the functionalities of scripts/build_offshore_shapes.py.
"""

import geopandas as gpd
import pytest
from shapely.geometry import box

from scripts._helpers import read_geo_boundaries
from scripts.build_offshore_shapes import build_offshore_shapes


@pytest.fixture
def geo_boundaries_parquet(tmp_path):
    shapes = gpd.GeoDataFrame(
        {
            "country_id": ["DEU", "DEU", "DEU", "FRA"],
            "shape_class": ["land", "maritime", "maritime", "maritime"],
            "parent": ["nuts", "eez", "eez", "eez"],
            "parent_id": ["DE1", "1", "2", "3"],
            "parent_name": ["Baden", "EEZ", "EEZ", "EEZ"],
        },
        geometry=[box(0, 0, 1, 1), box(1, 0, 2, 1), box(2, 0, 3, 1), box(0, 1, 1, 2)],
        crs="EPSG:4326",
    )
    path = tmp_path / "shapes.parquet"
    shapes.to_parquet(path)
    return path


def test_read_geo_boundaries(geo_boundaries_parquet):
    land = read_geo_boundaries(geo_boundaries_parquet, ["DE", "FR"], "land")
    assert land["country"].tolist() == ["DE"]
    maritime = read_geo_boundaries(geo_boundaries_parquet, ["DE"], "maritime")
    assert len(maritime) == 2


def test_build_offshore_shapes(geo_boundaries_parquet):
    maritime = read_geo_boundaries(geo_boundaries_parquet, ["DE", "FR"], "maritime")
    offshore = build_offshore_shapes(maritime)
    assert offshore.index.name == "name"
    assert sorted(offshore.index) == ["DE", "FR"]
    assert offshore.loc["DE", "geometry"].area == pytest.approx(2.0)
