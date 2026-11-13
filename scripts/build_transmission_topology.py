# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Build the candidate corridors for new H2 and CO2 pipelines and biomass transport between clustered buses.

The bus coordinates are triangulated in a metric projection (EPSG:3035). Each
edge gets its great-circle length and the fraction of its length that lies in
offshore regions. The selection below runs once per carrier with its own
settings; edges with an offshore length above the carrier's maximum are
excluded. The remaining edges of the Gabriel graph form the base set: an edge belongs to
it if the circle with the edge as its diameter contains no other bus [Gabriel
and Sokal, 1969]. As both ends lie on this circle, the test only needs the
distance from the edge midpoint to the nearest bus, found with a k-d tree
[Bentley, 1975], and only Delaunay edges need testing, as the Gabriel graph is
a subgraph of the Delaunay triangulation [Matula and Sokal, 1980]. Then the
shortest remaining Delaunay edges are added until each bus reaches the minimum
degree or has no more Delaunay neighbours. A minimum degree of zero keeps the
full Delaunay graph.

References
----------
- Gabriel and Sokal (1969), [A New Statistical Approach to Geographic Variation Analysis](https://doi.org/10.2307/2412323)
- Matula and Sokal (1980), [Properties of Gabriel Graphs Relevant to Geographic Variation Research and the Clustering of Points in the Plane](https://doi.org/10.1111/j.1538-4632.1980.tb00031.x)
- Bentley (1975), [Multidimensional binary search trees used for associative searching](https://doi.org/10.1145/361002.361007)
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
from pypsa.geo import haversine_pts
from scipy.spatial import Delaunay, KDTree
from shapely.geometry import LineString

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

DISTANCE_CRS = "EPSG:3035"


def delaunay_edges(coords: np.ndarray) -> np.ndarray:
    """
    Return the unique Delaunay edges as sorted index pairs ``(i, j)`` with ``i < j``.

    Up to three points are connected pairwise. Otherwise, the points are
    triangulated after a tiny deterministic displacement (Qhull option ``QJ``),
    so that collinear and coincident points are kept as vertices.
    """
    if len(coords) <= 3:
        return np.column_stack(np.triu_indices(len(coords), k=1))
    triangles = Delaunay(coords, qhull_options="QJ").simplices
    edges = np.vstack(
        [triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [0, 2]]]
    )
    return np.unique(np.sort(edges, axis=1), axis=0)


def is_gabriel(edges: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Flag the edges with no other point closer to their midpoint than half their length."""
    midpoints = coords[edges].mean(axis=1)
    radius = np.linalg.norm(coords[edges[:, 0]] - coords[edges[:, 1]], axis=1) / 2
    distance, _ = KDTree(coords).query(midpoints)
    return distance >= radius * (1 - 1e-9)


def build_delaunay_graph(
    buses: pd.DataFrame, offshore_shapes: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Build the Delaunay edges between buses.

    Parameters
    ----------
    buses : pd.DataFrame
        Buses with coordinates ``x`` and ``y`` in EPSG:4326.
    offshore_shapes : gpd.GeoDataFrame
        Offshore regions used for the offshore fraction of each edge.

    Returns
    -------
    gpd.GeoDataFrame
        One row per edge with columns name, bus0, bus1 (in lexicographic order),
        length (km), gabriel_edge, underwater_fraction and geometry.
    """
    lonlat = buses[["x", "y"]].to_numpy()
    points = gpd.GeoSeries(gpd.points_from_xy(*lonlat.T), crs="EPSG:4326")
    coords = points.to_crs(DISTANCE_CRS).get_coordinates().to_numpy()

    edges = delaunay_edges(coords)
    logger.info(f"Delaunay triangulation has {len(edges)} edges.")

    bus0, bus1 = np.sort(buses.index.to_numpy()[edges], axis=1).T
    graph = gpd.GeoDataFrame(
        {
            "name": bus0 + " -> " + bus1,
            "bus0": bus0,
            "bus1": bus1,
            "length": haversine_pts(lonlat[edges[:, 0]], lonlat[edges[:, 1]]),
            "gabriel_edge": is_gabriel(edges, coords),
        },
        geometry=[LineString(lonlat[e]) for e in edges],
        crs="EPSG:4326",
    )

    lines = graph.geometry.to_crs(DISTANCE_CRS)
    offshore = offshore_shapes.to_crs(DISTANCE_CRS).union_all()
    graph["underwater_fraction"] = (
        (lines.intersection(offshore).length / lines.length).fillna(0.0).round(2)
    )
    return graph


def select_corridors(
    graph: gpd.GeoDataFrame,
    buses: pd.Index,
    min_degree: int,
    max_offshore_distance: float,
) -> pd.Series:
    """
    Select the Gabriel edges and add the shortest other edges up to a minimum degree.

    Parameters
    ----------
    graph : gpd.GeoDataFrame
        Delaunay edges with columns bus0, bus1, length, underwater_fraction and
        gabriel_edge.
    buses : pd.Index
        All buses, used to report buses below the minimum degree.
    min_degree : int
        Minimum number of edges per bus. Values <= 0 select all edges.
    max_offshore_distance : float
        Maximum offshore length of an edge (km).

    Returns
    -------
    pd.Series
        Boolean mask of the selected edges.
    """
    allowed = graph["length"] * graph["underwater_fraction"] <= max_offshore_distance
    if min_degree <= 0:
        return allowed

    selected = (graph["gabriel_edge"] & allowed).to_numpy(copy=True)
    gabriel_ends = graph.loc[selected, ["bus0", "bus1"]].stack().value_counts()
    degree = pd.Series(0, index=buses).add(gabriel_ends, fill_value=0)

    for i in np.argsort(graph["length"].to_numpy(), kind="stable"):
        u, v = graph["bus0"].iat[i], graph["bus1"].iat[i]
        if (
            allowed.iat[i]
            and not selected[i]
            and min(degree[u], degree[v]) < min_degree
        ):
            selected[i] = True
            degree[[u, v]] += 1

    if not (unmet := degree.index[degree < min_degree]).empty:
        logger.warning(
            f"{len(unmet)} of {len(degree)} buses have fewer than {min_degree} "
            f"candidate corridors: {', '.join(unmet)}"
        )
    return pd.Series(selected, index=graph.index)


def build_transmission_topology(
    n: pypsa.Network,
    offshore_shapes: gpd.GeoDataFrame,
    corridors: dict[str, dict],
) -> gpd.GeoDataFrame:
    """
    Build the candidate corridors; see the module docstring for the method.

    Parameters
    ----------
    n : pypsa.Network
        Network with bus coordinates ``x`` and ``y`` in EPSG:4326.
    offshore_shapes : gpd.GeoDataFrame
        Offshore regions used for the offshore fraction of each edge.
    corridors : dict[str, dict]
        Per carrier, the ``min_degree`` and ``max_offshore_distance``.

    Returns
    -------
    gpd.GeoDataFrame
        Delaunay edges selected for at least one carrier, with one boolean
        column per carrier.
    """
    buses = n.buses.query("carrier == 'AC'")
    graph = build_delaunay_graph(buses, offshore_shapes)
    for carrier, settings in corridors.items():
        logger.info(f"Selecting candidate corridors for {carrier}.")
        graph[carrier] = select_corridors(graph, buses.index, **settings)
    return graph[graph[list(corridors)].any(axis=1)]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_transmission_topology")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes)
    candidates = build_transmission_topology(
        n, offshore_shapes, snakemake.params.transport_corridors
    )
    candidates.to_file(snakemake.output.candidates)
