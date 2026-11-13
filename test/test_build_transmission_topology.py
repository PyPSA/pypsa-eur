# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""Tests for the candidate transport corridors in scripts/build_transmission_topology.py."""

import geopandas as gpd
import pypsa
import pytest
from shapely.geometry import box

from scripts.build_transmission_topology import build_transmission_topology

OFFSHORE = gpd.GeoDataFrame(geometry=[box(30, 30, 31, 31)], crs="EPSG:4326")


def network(coords: list[tuple[float, float]]) -> pypsa.Network:
    n = pypsa.Network()
    for i, (x, y) in enumerate(coords):
        n.add("Bus", f"b{i}", x=x, y=y, carrier="AC")
    n.add("Bus", "b0 H2", x=0.0, y=50.0, carrier="H2")
    return n


@pytest.mark.parametrize(
    "coords, min_degree, expected",
    [
        ([(10, 50)], 1, set()),
        ([(10, 50), (11, 50)], 1, {("b0", "b1")}),
        (
            [(10, 50), (11, 50), (12, 50), (13, 50)],
            1,
            {("b0", "b1"), ("b1", "b2"), ("b2", "b3")},
        ),
        (
            [(10, 50), (12, 50), (11, 51)],
            1,
            {("b0", "b1"), ("b0", "b2"), ("b1", "b2")},
        ),
        (
            [(10, 50), (11, 50), (11, 51), (10, 51)],
            1,
            {("b0", "b1"), ("b1", "b2"), ("b2", "b3"), ("b0", "b3")},
        ),
        (
            [(10, 50), (10, 50), (11, 50)],
            1,
            {("b0", "b1"), ("b0", "b2"), ("b1", "b2")},
        ),
        ([(10, 50), (12, 50), (11, 50.1)], 1, {("b0", "b2"), ("b1", "b2")}),
        (
            [(10, 50), (12, 50), (11, 50.1)],
            2,
            {("b0", "b1"), ("b0", "b2"), ("b1", "b2")},
        ),
        (
            [(10, 50), (12, 50), (11, 50.1)],
            0,
            {("b0", "b1"), ("b0", "b2"), ("b1", "b2")},
        ),
    ],
    ids=[
        "point",
        "pair",
        "collinear",
        "triangle",
        "square",
        "coincident",
        "obtuse-gabriel",
        "obtuse-min-degree",
        "obtuse-delaunay",
    ],
)
def test_candidate_corridors(coords, min_degree, expected):
    """Small and degenerate layouts give the expected corridors between AC buses."""
    settings = {"min_degree": min_degree, "max_offshore_distance": float("inf")}
    candidates = build_transmission_topology(
        network(coords), OFFSHORE, {"H2": settings}
    )
    assert set(zip(candidates.bus0, candidates.bus1)) == expected


def test_carrier_specific_corridors():
    """Each carrier gets its own selection; the union of all is returned."""
    inf = float("inf")
    corridors = {
        "H2": {"min_degree": 1, "max_offshore_distance": inf},
        "CO2": {"min_degree": 0, "max_offshore_distance": inf},
    }
    candidates = build_transmission_topology(
        network([(10, 50), (12, 50), (11, 50.1)]), OFFSHORE, corridors
    ).set_index(["bus0", "bus1"])
    assert candidates["H2"].to_dict() == {
        ("b0", "b1"): False,
        ("b0", "b2"): True,
        ("b1", "b2"): True,
    }
    assert candidates["CO2"].all()
