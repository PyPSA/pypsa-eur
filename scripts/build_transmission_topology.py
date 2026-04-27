# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Filter the Delaunay graph to transmission topology candidates.

Description
-----------
Reads the full Delaunay edge table produced by build_transmission_delaunay_graph
and applies Gabriel graph filtering and min-degree backfilling to produce a
candidate edge set for one (min_degree, max_offshore_haversine_distance) pair.

Outputs
-------
GeoJSON file in WGS84 (EPSG:4326) format:

``candidates``
    Selected candidate corridor table (subset of Delaunay edges) with columns:

    - ``name``
    - ``bus0``
    - ``bus1``
    - ``length``
    - ``gabriel_edge``
    - ``underwater_fraction``
    - ``geometry``
"""

import logging

import geopandas as gpd
import numpy as np

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

COLS_EDGES = [
    "name",
    "bus0",
    "bus1",
    "length",
    "gabriel_edge",
    "geometry",
    "underwater_fraction",
]


def enforce_min_degree(
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
    delaunay_graph : gpd.GeoDataFrame
        Full Delaunay triangulation edge table.
        Must contain integer columns ``source``, ``target``, and ``length``,
        with canonical ordering ``source < target`` and a boolean
        ``gabriel_edge`` column.
    node_count : int
        Total number of nodes in the graph. Used to size degree arrays.
    min_degree : int
        Minimum number of edges required per node. Values <= 0 disable
        Gabriel filtering and return the full ``delaunay_graph``.

    Returns
    -------
    gpd.GeoDataFrame
        Edge GeoDataFrame with the same schema as ``delaunay_graph``.
        For ``min_degree > 0``, starts from Gabriel edges and extends with any
        backfill edges needed to meet the degree constraint.
        Nodes whose maximum possible Delaunay degree is below ``min_degree``
        are handled gracefully with a warning.

    Warns
    -----
    Logs a warning if the minimum degree cannot be fully satisfied, either
    because a node has too few Delaunay neighbours (unreachable) or because
    the greedy backfill did not converge (unmet after backfill).
    """
    if min_degree <= 0:
        logger.info(
            "Gabriel filtering disabled via min_degree=%s; returning full Delaunay graph.",
            min_degree,
        )
        return delaunay_graph.copy()

    import pandas as pd

    result = delaunay_graph[delaunay_graph["gabriel_edge"]].copy()
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


def filter_by_max_offshore_haversine_distance(
    delaunay_graph: gpd.GeoDataFrame,
    max_offdist: float = float("inf"),
) -> gpd.GeoDataFrame:
    """
    Filter Delaunay edges by maximum offshore haversine distance.

    Parameters
    ----------
    delaunay_graph : gpd.GeoDataFrame
        Full Delaunay edge table.
    max_offdist : float, default inf
        Maximum allowed offshore edge length in km, computed as
        ``underwater_fraction * length``. ``inf`` disables this filter.

    Returns
    -------
    gpd.GeoDataFrame
        Delaunay edge table after applying offshore-distance filtering.
    """
    if np.isfinite(max_offdist):
        offshore_length = (
            delaunay_graph["length"] * delaunay_graph["underwater_fraction"]
        )
        return delaunay_graph[offshore_length <= max_offdist].copy()
    return delaunay_graph.copy()


if __name__ == "__main__":
    is_mock_run = "snakemake" not in globals()

    if is_mock_run:
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_transmission_topology",
            clusters="50",
            min_degree=1,
            max_offdist="inf",
            configfiles=["config/config.200.yaml"],
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    min_degree = snakemake.params["min_degree"]
    max_offdist = snakemake.params["max_offdist"]

    delaunay_graph = gpd.read_file(snakemake.input.delaunay_graph)

    source_arr = delaunay_graph["source"].to_numpy(dtype=int)
    target_arr = delaunay_graph["target"].to_numpy(dtype=int)
    node_count = int(max(source_arr.max(), target_arr.max())) + 1

    logger.info(
        "Filtering topology: min_degree=%s, max_offdist=%s.",
        min_degree,
        max_offdist,
    )

    candidate_pool = filter_by_max_offshore_haversine_distance(
        delaunay_graph=delaunay_graph,
        max_offdist=max_offdist,
    )

    selected_edges = enforce_min_degree(
        delaunay_graph=candidate_pool,
        node_count=node_count,
        min_degree=min_degree,
    )[COLS_EDGES]

    selected_edges.to_file(snakemake.output.candidates)
