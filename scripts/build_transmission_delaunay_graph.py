# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Build the Delaunay triangulation graph for transmission corridor candidates.

Description
-----------
Creates the full Delaunay edge table from clustered network bus coordinates,
enriched with Gabriel edge flags and offshore underwater fractions.

The output GeoDataFrame is consumed by build_transmission_topology to filter
candidate edges per carrier configuration.

Outputs
-------
GeoJSON file in WGS84 (EPSG:4326) format:

``delaunay_graph``
    Full Delaunay edge table with columns:

    - ``name``: Canonical undirected edge identifier ``"bus0 -> bus1"``.
    - ``bus0``: Canonically ordered first bus id (lexicographic order).
    - ``bus1``: Canonically ordered second bus id.
    - ``source``: Integer index of the first bus.
    - ``target``: Integer index of the second bus.
    - ``length``: Great-circle edge length in km.
    - ``gabriel_edge``: ``True`` if the edge satisfies the Gabriel empty
        circle criterion.
    - ``underwater_fraction``: Fraction of edge length over offshore regions.
    - ``geometry``: LineString geometry in ``EPSG:4326``.
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
from numpy.typing import NDArray
from pypsa.geo import haversine_pts
from scipy.spatial import Delaunay, QhullError
from shapely.geometry import LineString

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

GEO_CRS = "EPSG:4326"
DISTANCE_CRS = "EPSG:3035"


def delaunay_edges(coords_meter: NDArray[np.float64]) -> list[tuple[int, int]]:
    """
    Build sorted unique undirected edges from Delaunay simplices.

    Parameters
    ----------
    coords_meter : NDArray[np.float64]
        Node coordinates in a metric CRS with shape ``(n_nodes, 2)``.

    Returns
    -------
    list[tuple[int, int]]
        Sorted list of unique undirected edges represented as index pairs
        ``(u, v)`` with ``u < v``.

    Raises
    ------
    ValueError
        If fewer than three points are provided.
    """
    if len(coords_meter) < 3:
        raise ValueError("Need at least 3 points to compute Delaunay triangulation.")

    try:
        tri = Delaunay(coords_meter)
    except QhullError:
        # Robust fallback for nearly degenerate point sets.
        tri = Delaunay(coords_meter, qhull_options="QJ")

    edges: set[tuple[int, int]] = set()
    for simplex in tri.simplices:
        for i, j in ((0, 1), (1, 2), (0, 2)):
            edges.add(tuple(sorted((int(simplex[i]), int(simplex[j])))))
    return sorted(edges)


def pairwise_sq_dists(coords: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Compute the squared Euclidean distance matrix.

    Parameters
    ----------
    coords : NDArray[np.float64]
        Node coordinates in Euclidean space with shape ``(n_nodes, n_dims)``.

    Returns
    -------
    NDArray[np.float64]
        Dense matrix ``D`` with shape ``(n_nodes, n_nodes)`` where
        ``D[i, j] = ||coords[i] - coords[j]||^2``.
    """
    norms = np.sum(coords * coords, axis=1)
    sq = norms[:, None] + norms[None, :] - 2.0 * (coords @ coords.T)
    return np.maximum(sq, 0.0)


def classify_gabriel_edges(
    edges: NDArray[np.int64],
    coords_meter: NDArray[np.float64],
    eps: float = 1e-9,
    chunk_size: int = 1024,
) -> NDArray[np.bool_]:
    """
    Classify candidate edges as Gabriel or non-Gabriel.

    Parameters
    ----------
    edges : NDArray[np.int64]
        Candidate undirected edges as node index pairs with shape ``(n_edges, 2)``.
    coords_meter : NDArray[np.float64]
        Node coordinates in a metric CRS with shape ``(n_nodes, 2)``.
    eps : float, default 1e-9
        Numerical tolerance used in the Gabriel emptiness test.
    chunk_size : int, default 1024
        Batch size used to limit temporary memory during vectorized checks.

    Returns
    -------
    NDArray[np.bool_]
        Boolean mask of shape ``(n_edges,)`` where ``True`` marks Gabriel edges.
    """
    sq_dists = pairwise_sq_dists(coords_meter)
    edge_count = edges.shape[0]
    is_gabriel = np.ones(edge_count, dtype=bool)

    for start in range(0, edge_count, chunk_size):
        end = min(start + chunk_size, edge_count)
        uv = edges[start:end]
        u = uv[:, 0]
        v = uv[:, 1]
        dij2 = sq_dists[u, v]

        # A point k invalidates edge (u, v) when dik^2 + djk^2 < dij^2.
        inside = (sq_dists[u, :] + sq_dists[v, :]) < (dij2[:, None] - eps)
        row = np.arange(end - start)
        inside[row, u] = False
        inside[row, v] = False
        is_gabriel[start:end] = ~inside.any(axis=1)

    return is_gabriel


def get_bus_coordinates(
    n: pypsa.Network,
) -> tuple[pd.Index, NDArray[np.float64], NDArray[np.float64]]:
    """
    Extract bus ids and coordinates in geographic and metric CRS.

    Parameters
    ----------
    n : pypsa.Network
        Input network containing bus coordinate columns ``x`` and ``y``.

    Returns
    -------
    tuple[pd.Index, NDArray[np.float64], NDArray[np.float64]]
        Tuple ``(bus_ids, coords_geo, coords_meter)`` where ``coords_geo`` is in
        ``EPSG:4326`` and ``coords_meter`` is in ``EPSG:3035``.

    Raises
    ------
    ValueError
        If fewer than two buses with valid coordinates are available.
    """
    buses = n.buses.copy()
    buses = buses.dropna(subset=["x", "y"])

    if len(buses) < 2:
        raise ValueError("Need at least 2 buses with coordinates to build corridors.")

    points_geo = gpd.points_from_xy(buses["x"], buses["y"], crs=GEO_CRS)
    points_meter = gpd.GeoSeries(points_geo, crs=GEO_CRS).to_crs(DISTANCE_CRS)

    coords_geo = np.column_stack(
        (np.asarray(points_geo.x), np.asarray(points_geo.y))
    ).astype(float)
    coords_meter = np.column_stack(
        (np.asarray(points_meter.x), np.asarray(points_meter.y))
    ).astype(float)

    return buses.index, coords_geo, coords_meter


def delaunay_triangulation(
    bus_ids: pd.Index,
    coords_geo: NDArray[np.float64],
    coords_meter: NDArray[np.float64],
) -> gpd.GeoDataFrame:
    """
    Build the full Delaunay edge table with geometry and metadata.

    Parameters
    ----------
    bus_ids : pd.Index
        Bus identifiers aligned with the coordinate arrays.
    coords_geo : NDArray[np.float64]
        Geographic coordinates (longitude, latitude) in ``EPSG:4326``.
    coords_meter : NDArray[np.float64]
        Projected coordinates in ``EPSG:3035`` used for metric calculations.

    Returns
    -------
    gpd.GeoDataFrame
        Delaunay edge table including integer node endpoints, normalized bus
        names, canonical edge ``name``, ``length`` (haversine distance, km),
        ``gabriel_edge`` flag, and LineString geometry.
    """

    edges = delaunay_edges(coords_meter)
    edges_arr = np.array(edges, dtype=np.int64)
    u = edges_arr[:, 0]
    v = edges_arr[:, 1]

    logger.info("Compute Delaunay triangulation with %s edges.", len(edges))
    lengths = haversine_pts(coords_geo[u], coords_geo[v])
    gabriel_flags = classify_gabriel_edges(edges_arr, coords_meter)
    edge_geoms = [LineString([coords_geo[i], coords_geo[j]]) for i, j in edges]

    bus0 = bus_ids.take(u).astype(str).to_numpy()
    bus1 = bus_ids.take(v).astype(str).to_numpy()

    # Keep undirected edge naming compatible with topology naming in
    # prepare_sector_network.create_network_topology.
    swap_mask = bus0 > bus1
    bus0_norm = np.where(swap_mask, bus1, bus0)
    bus1_norm = np.where(swap_mask, bus0, bus1)
    connector = " -> "
    edge_names = np.char.add(
        np.char.add(bus0_norm, connector),
        bus1_norm,
    )

    edges_gdf = gpd.GeoDataFrame(
        {
            "source": u,
            "target": v,
            "bus0": bus0_norm,
            "bus1": bus1_norm,
            "name": edge_names,
            "length": lengths,
            "gabriel_edge": gabriel_flags,
            "geometry": edge_geoms,
        },
        geometry="geometry",
        crs=GEO_CRS,
    )
    return edges_gdf


def add_underwater_fraction(
    edges: gpd.GeoDataFrame,
    offshore_shapes: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Compute the fraction of each edge's length that lies within offshore regions.

    Parameters
    ----------
    edges : gpd.GeoDataFrame
        Edge table with LineString geometry in ``EPSG:4326``.
    offshore_shapes : gpd.GeoDataFrame
        Offshore region polygons used to determine submarine sections.

    Returns
    -------
    gpd.GeoDataFrame
        ``edges`` with an added ``underwater_fraction`` column in ``[0, 1]``.
    """
    offshore_union = offshore_shapes.to_crs(DISTANCE_CRS).union_all()
    geoms = edges.geometry.to_crs(DISTANCE_CRS)
    edges["underwater_fraction"] = (
        (geoms.intersection(offshore_union).length / geoms.length).fillna(0.0).round(2)
    )
    return edges


if __name__ == "__main__":
    is_mock_run = "snakemake" not in globals()

    if is_mock_run:
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_transmission_delaunay_graph",
            clusters="50",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes)

    bus_ids, coords_geo, coords_meter = get_bus_coordinates(n)
    delaunay_graph = delaunay_triangulation(
        bus_ids=bus_ids,
        coords_geo=coords_geo,
        coords_meter=coords_meter,
    )

    add_underwater_fraction(delaunay_graph, offshore_shapes)

    delaunay_graph.to_file(snakemake.output.delaunay_graph)
