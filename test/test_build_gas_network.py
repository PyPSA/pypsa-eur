# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>>
#
# SPDX-License-Identifier: MIT

"""
Tests the functionalities of scripts/build_gas_network.py.
"""

import geopandas as gpd
import pytest
from shapely.geometry import LineString

from scripts.build_gas_network import diameter_to_capacity, prepare_dataset

MCM_PER_DAY_TO_MW = 437.5


@pytest.mark.parametrize(
    "p_nom,length,p_nom_expected,length_expected",
    [
        # plausible capacity and length are kept
        (20000.0, 120.0, 20000.0, 120.0),
        # implausible capacity is replaced by the diameter-based capacity
        (300.0, 120.0, diameter_to_capacity(1000.0), 120.0),
        # implausible length is replaced by the scaled haversine distance
        (20000.0, 900.0, 20000.0, pytest.approx(107.0, rel=0.01)),
    ],
)
def test_prepare_dataset_corrections(p_nom, length, p_nom_expected, length_expected):
    """
    Verify that prepare_dataset only corrects implausible capacities and lengths.
    """
    df = gpd.GeoDataFrame(
        {
            "name": ["pipe"],
            "diameter_mm": [1000.0],
            "diameter_mm_method": ["raw"],
            "is_H_gas": [1],
            "is_bothDirection": [0],
            "length_km": [length],
            "max_cap_M_m3_per_d": [p_nom / MCM_PER_DAY_TO_MW],
            "max_pressure_bar": [80.0],
            "start_year": [2000],
        },
        geometry=[LineString([(10.0, 50.0), (11.0, 50.0)])],
        crs=4326,
    )

    df = prepare_dataset(df)

    assert df.at[0, "p_nom"] == pytest.approx(p_nom_expected)
    assert df.at[0, "length"] == length_expected
