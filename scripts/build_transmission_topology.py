# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Build transmission corridor topology/candidates from clustered network
bus coordinates using Delaunay triangulation and optional Gabriel graph
filtering.

Description
-----------
Creates the topology for potential investment corridors between
clustered buses, i.e. for H2 and CO2 pipeline candidates.

The script computes the Delaunay triangulation of the bus coordinates
and optionally filters edges to retain only Gabriel edges, which are
more likely to represent direct connections. A minimum degree
constraint can be enforced to ensure network connectivity and to avoid
stubs.

Outputs
-------
Two GeoJSON files in WGS84 (EPSG:4326) format are written per carrier and cluster configuration:

1. ``all_edges``
     Full Delaunay edge table (one row per Delaunay edge) with columns:

     - ``name``: Canonical undirected edge identifier ``"bus0 -> bus1"``.
     - ``bus0``: Canonically ordered first bus id (lexicographic order).
     - ``bus1``: Canonically ordered second bus id.
     - ``length``: Great-circle edge length in km multiplied by
         ``length_factor``.
     - ``gabriel_edge``: ``True`` if the edge satisfies the Gabriel empty
         circle criterion.
     - ``selected_edge``: ``True`` if the edge is part of the final
         candidate set after filtering/backfilling.
     - ``geometry``: LineString geometry in ``EPSG:4326``.

2. ``candidates``
     Selected candidate corridor table (subset of Delaunay edges) with
     columns:

     - ``name``
     - ``bus0``
     - ``bus1``
     - ``length``
     - ``gabriel_edge``
     - ``geometry``

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
COLS_EDGES = [
    "name",
    "bus0",
    "bus1",
    "length",
    "gabriel_edge",
    "geometry",
]
BUS_SUFFIX = {
    "carbon_dioxide": "co2 stored",
    "hydrogen": "H2",
}
LINK_PREFIX = {
    "carbon_dioxide": "CO2 pipeline",
    "hydrogen": "H2 pipeline",
}


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
    length_factor: float = 1.25,
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
    logger.info("Calculating edge lengths using length factor %.2f.", length_factor)
    lengths = haversine_pts(coords_geo[u], coords_geo[v]) * length_factor
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


def enforce_min_degree(
    selected_edges: gpd.GeoDataFrame,
    delaunay_graph: gpd.GeoDataFrame,
    node_count: int,
    min_degree: int,
) -> gpd.GeoDataFrame:
    """
    Add shortest Delaunay edges to satisfy a minimum node degree.

    Iterates over all Delaunay edges in ascending order of length and adds
    edges to the selected set until every node reaches the requested minimum
    degree, or until no further Delaunay edges are available for that node.

    Parameters
    ----------
    selected_edges : gpd.GeoDataFrame
        Currently selected candidate edges, e.g. after Gabriel filtering.
        Must contain integer columns ``source``, ``target``, and ``length``.
    delaunay_graph : gpd.GeoDataFrame
        Full Delaunay triangulation edges used as the backfill pool.
        Must contain integer columns ``source``, ``target``, and ``length``,
        with canonical ordering ``source < target``.
    node_count : int
        Total number of nodes in the graph. Used to size degree arrays.
    min_degree : int
        Minimum number of edges required per node. Values <= 0 disable the
        constraint and return ``selected_edges`` unchanged.

    Returns
    -------
    gpd.GeoDataFrame
        Edge GeoDataFrame with the same schema as ``selected_edges``,
        extended with any backfill edges needed to meet the degree constraint.
        Nodes whose maximum possible Delaunay degree is below ``min_degree``
        are handled gracefully with a warning.

    Warns
    -----
    Logs a warning if the minimum degree cannot be fully satisfied, either
    because a node has too few Delaunay neighbours (unreachable) or because
    the greedy backfill did not converge (unmet after backfill).
    """
    if min_degree <= 0:
        logger.info("Min-degree constraint disabled (min_degree=%s).", min_degree)
        return selected_edges

    result = selected_edges.copy()
    all_sorted = delaunay_graph.sort_values("length").reset_index(drop=True)
    source_arr = all_sorted["source"].to_numpy(dtype=int)
    target_arr = all_sorted["target"].to_numpy(dtype=int)

    max_possible_degree = np.bincount(
        np.concatenate([source_arr, target_arr]), minlength=node_count
    )

    sel_source = result["source"].to_numpy(dtype=int)
    sel_target = result["target"].to_numpy(dtype=int)
    edge_set: set[tuple[int, int]] = set(zip(sel_source, sel_target))
    degrees = np.bincount(
        np.concatenate([sel_source, sel_target]), minlength=node_count
    )

    deficits = np.maximum(0, min_degree - degrees)
    remaining_deficits = int(deficits.sum())
    additions_idx: list[int] = []

    for idx, (u, v) in enumerate(zip(source_arr, target_arr)):
        if remaining_deficits == 0:
            break
        if (u, v) in edge_set:
            continue
        if deficits[u] <= 0 and deficits[v] <= 0:
            continue
        edge_set.add((u, v))
        additions_idx.append(idx)
        if deficits[u] > 0:
            deficits[u] -= 1
            remaining_deficits -= 1
        if deficits[v] > 0:
            deficits[v] -= 1
            remaining_deficits -= 1

    if additions_idx:
        to_add = all_sorted.iloc[additions_idx]
        result = gpd.GeoDataFrame(
            pd.concat([result, to_add], ignore_index=True),
            geometry="geometry",
            crs=delaunay_graph.crs,
        )

    result_source = result["source"].to_numpy(dtype=int)
    result_target = result["target"].to_numpy(dtype=int)
    achieved_degree = np.bincount(
        np.concatenate([result_source, result_target]), minlength=node_count
    )
    target_degree = np.minimum(min_degree, max_possible_degree)
    unreachable_mask = max_possible_degree < min_degree
    unmet_mask = achieved_degree < target_degree

    if np.any(unreachable_mask) or np.any(unmet_mask):
        logger.warning(
            "Min-degree not fully met: requested=%s, unreachable=%s/%s, unmet_after_backfill=%s/%s.",
            min_degree,
            int(unreachable_mask.sum()),
            node_count,
            int(unmet_mask.sum()),
            node_count,
        )

    return result


def prepare_candidate_edges(
    delaunay_graph: gpd.GeoDataFrame,
    gabriel_filter_enabled: bool,
    node_count: int,
    min_degree: int,
) -> gpd.GeoDataFrame:
    """
    Build selected candidate edges from a Delaunay edge table.

    Parameters
    ----------
    delaunay_graph : gpd.GeoDataFrame
        Full Delaunay edge table.
    gabriel_filter_enabled : bool
        Whether to restrict the initial selection to Gabriel edges.
    node_count : int
        Number of nodes used for degree accounting.
    min_degree : int
        Minimum per-node degree target used for optional backfilling.

    Returns
    -------
    gpd.GeoDataFrame
        Selected edge table after optional Gabriel filtering and min-degree
        backfilling.
    """
    if not gabriel_filter_enabled:
        return delaunay_graph.copy()

    selected_edges = delaunay_graph[delaunay_graph["gabriel_edge"]].copy()
    selected = enforce_min_degree(
        selected_edges=selected_edges,
        delaunay_graph=delaunay_graph,
        node_count=node_count,
        min_degree=min_degree,
    )

    return selected


def mark_selected_edges(
    delaunay_graph: gpd.GeoDataFrame,
    selected_edges: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Annotate all Delaunay edges with selection membership.

    Parameters
    ----------
    delaunay_graph : gpd.GeoDataFrame
        Full Delaunay edge table to annotate.
    selected_edges : gpd.GeoDataFrame
        Selected candidate edges containing the ``name`` column.

    Returns
    -------
    gpd.GeoDataFrame
        ``delaunay_graph`` with an added boolean ``selected_edge`` column.
    """
    selected_names = set(selected_edges["name"])
    delaunay_graph["selected_edge"] = delaunay_graph["name"].isin(selected_names)
    return delaunay_graph


if __name__ == "__main__":
    is_mock_run = "snakemake" not in globals()

    if is_mock_run:
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_transmission_topology",
            carrier="hydrogen",
            clusters="200",
            run="nodes200",
            configfiles=["config/config.200.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    carrier_enabled = bool(snakemake.params["carrier_enabled"])
    carrier = snakemake.wildcards.carrier
    gabriel_filter_enabled = bool(snakemake.params["gabriel_filter_enabled"])
    min_degree = int(snakemake.params["min_degree"])
    length_factor = float(snakemake.params["length_factor"])

    if not carrier_enabled:
        raise ValueError(
            "Transmission carrier is disabled in config for this wildcard."
        )

    n = pypsa.Network(snakemake.input.network)

    # Delaunay triangulation and Gabriel filtering
    bus_ids, coords_geo, coords_meter = get_bus_coordinates(n)
    delaunay_graph = delaunay_triangulation(
        bus_ids=bus_ids,
        coords_geo=coords_geo,
        coords_meter=coords_meter,
        length_factor=length_factor,
    )
    selected_edges = prepare_candidate_edges(
        delaunay_graph=delaunay_graph,
        gabriel_filter_enabled=gabriel_filter_enabled,
        node_count=len(coords_geo),
        min_degree=min_degree,
    )
    selected_edges = selected_edges[COLS_EDGES]
    delaunay_graph = mark_selected_edges(delaunay_graph, selected_edges)

    # Rename buses and link names
    delaunay_graph[["bus0", "bus1"]] = (
        delaunay_graph[["bus0", "bus1"]] + " " + BUS_SUFFIX[carrier]
    )
    delaunay_graph["name"] = LINK_PREFIX[carrier] + " " + delaunay_graph["name"]
    selected_edges[["bus0", "bus1"]] = (
        selected_edges[["bus0", "bus1"]] + " " + BUS_SUFFIX[carrier]
    )
    selected_edges["name"] = LINK_PREFIX[carrier] + " " + selected_edges["name"]

    # Export
    delaunay_graph.to_file(snakemake.output.all_edges)
    selected_edges.to_file(snakemake.output.candidates)
