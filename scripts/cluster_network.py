# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Creates networks clustered to ``{cluster}`` number of zones with aggregated
buses and transmission corridors.

Outputs
-------

- ``resources/regions_onshore_base_s_{clusters}.geojson``:

    .. image:: img/regions_onshore_base_s_X.png
        :scale: 33 %

- ``resources/regions_offshore_base_s_{clusters}.geojson``:

    .. image:: img/regions_offshore_base_s_X.png
        :scale: 33 %

- ``resources/busmap_base_s_{clusters}.csv``: Mapping of buses from ``networks/base.nc`` to ``networks/base_s_{clusters}.nc``;
- ``resources/linemap_base_s_{clusters}.csv``: Mapping of lines from ``networks/base.nc`` to ``networks/base_s_{clusters}.nc``;
- ``networks/base_s_{clusters}.nc``:

    .. image:: img/base_s_X.png
        :scale: 40  %

Description
-----------

.. note::

    **Is it possible to run the model without the** ``simplify_network`` **rule?**

        No, the network clustering methods in the PyPSA module
        `pypsa.clustering.spatial <https://github.com/PyPSA/PyPSA/blob/master/pypsa/clustering/spatial.py>`_
        do not work reliably with multiple voltage levels and transformers.

Exemplary unsolved network clustered to 512 nodes:

.. image:: img/base_s_512.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 256 nodes:

.. image:: img/base_s_256.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 128 nodes:

.. image:: img/base_s_128.png
    :scale: 40  %
    :align: center

Exemplary unsolved network clustered to 37 nodes:

.. image:: img/base_s_37.png
    :scale: 40  %
    :align: center
"""

import logging
import warnings
from functools import reduce

import geopandas as gpd
import linopy
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import configure_logging, set_scenario_config
from packaging.version import Version, parse
from pypsa.clustering.spatial import (
    busmap_by_greedy_modularity,
    busmap_by_hac,
    busmap_by_kmeans,
    get_clustering_from_busmap,
)
from scipy.sparse.csgraph import connected_components

PD_GE_2_2 = parse(pd.__version__) >= Version("2.2")

warnings.filterwarnings(action="ignore", category=UserWarning)
idx = pd.IndexSlice
logger = logging.getLogger(__name__)


def normed(x):
    return (x / x.sum()).fillna(0.0)


def weighting_for_country(df: pd.DataFrame, weights: pd.Series) -> pd.Series:
    w = normed(weights.reindex(df.index, fill_value=0))
    return (w * (100 / w.max())).clip(lower=1).astype(int)


def get_feature_data_for_hac(fn: str) -> pd.DataFrame:
    ds = xr.open_dataset(fn)
    feature_data = (
        pd.concat([ds[var].to_pandas() for var in ds.data_vars], axis=0).fillna(0.0).T
    )
    feature_data.columns = feature_data.columns.astype(str)
    return feature_data


def fix_country_assignment_for_hac(n: pypsa.Network) -> None:
    # overwrite country of nodes that are disconnected from their country-topology
    for country in n.buses.country.unique():
        m = n[n.buses.country == country].copy()

        _, labels = connected_components(m.adjacency_matrix(), directed=False)

        component = pd.Series(labels, index=m.buses.index)
        component_sizes = component.value_counts()

        if len(component_sizes) > 1:
            disconnected_bus = component[component == component_sizes.index[-1]].index[
                0
            ]

            neighbor_bus = n.lines.query(
                "bus0 == @disconnected_bus or bus1 == @disconnected_bus"
            ).iloc[0][["bus0", "bus1"]]
            new_country = list(set(n.buses.loc[neighbor_bus].country) - {country})[0]

            logger.info(
                f"overwriting country `{country}` of bus `{disconnected_bus}` "
                f"to new country `{new_country}`, because it is disconnected "
                "from its initial inter-country transmission grid."
            )
            n.buses.at[disconnected_bus, "country"] = new_country


def distribute_n_clusters_to_countries(
    n: pypsa.Network,
    n_clusters: int,
    cluster_weights: pd.Series,
    focus_weights: dict | None = None,
    solver_name: str = "scip",
) -> pd.Series:
    """
    Determine the number of clusters per country.
    """
    L = (
        cluster_weights.groupby([n.buses.country, n.buses.sub_network])
        .sum()
        .pipe(normed)
    )

    N = n.buses.groupby(["country", "sub_network"]).size()[L.index]

    assert n_clusters >= len(N) and n_clusters <= N.sum(), (
        f"Number of clusters must be {len(N)} <= n_clusters <= {N.sum()} for this selection of countries."
    )

    if isinstance(focus_weights, dict):
        total_focus = sum(list(focus_weights.values()))

        assert total_focus <= 1.0, (
            "The sum of focus weights must be less than or equal to 1."
        )

        for country, weight in focus_weights.items():
            L[country] = weight / len(L[country])

        remainder = [
            c not in focus_weights.keys() for c in L.index.get_level_values("country")
        ]
        L[remainder] = L.loc[remainder].pipe(normed) * (1 - total_focus)

        logger.warning("Using custom focus weights for determining number of clusters.")

    assert np.isclose(L.sum(), 1.0, rtol=1e-3), (
        f"Country weights L must sum up to 1.0 when distributing clusters. Is {L.sum()}."
    )

    m = linopy.Model()
    clusters = m.add_variables(
        lower=1, upper=N, coords=[L.index], name="n", integer=True
    )
    m.add_constraints(clusters.sum() == n_clusters, name="tot")
    # leave out constant in objective (L * n_clusters) ** 2
    m.objective = (clusters * clusters - 2 * clusters * L * n_clusters).sum()
    if solver_name == "gurobi":
        logging.getLogger("gurobipy").propagate = False
    elif solver_name not in ["scip", "cplex", "xpress", "copt", "mosek"]:
        logger.info(
            f"The configured solver `{solver_name}` does not support quadratic objectives. Falling back to `scip`."
        )
        solver_name = "scip"
    m.solve(solver_name=solver_name)
    return m.solution["n"].to_series().astype(int)


def busmap_for_n_clusters(
    n: pypsa.Network,
    n_clusters_c: pd.Series,
    cluster_weights: pd.Series,
    algorithm: str = "kmeans",
    features: pd.DataFrame | None = None,
    **algorithm_kwds,
) -> pd.Series:
    if algorithm == "hac" and features is None:
        raise ValueError("For HAC clustering, features must be provided.")

    if algorithm == "kmeans":
        algorithm_kwds.setdefault("n_init", 1000)
        algorithm_kwds.setdefault("max_iter", 30000)
        algorithm_kwds.setdefault("tol", 1e-6)
        algorithm_kwds.setdefault("random_state", 0)

    def busmap_for_country(x):
        prefix = x.name[0] + x.name[1] + " "
        logger.debug(
            f"Determining busmap for country {prefix[:-1]} "
            f"from {len(x)} buses to {n_clusters_c[x.name]}."
        )
        if len(x) == 1:
            return pd.Series(prefix + "0", index=x.index)
        weight = weighting_for_country(x, cluster_weights)

        if algorithm == "kmeans":
            return prefix + busmap_by_kmeans(
                n, weight, n_clusters_c[x.name], buses_i=x.index, **algorithm_kwds
            )
        elif algorithm == "hac":
            return prefix + busmap_by_hac(
                n,
                n_clusters_c[x.name],
                buses_i=x.index,
                feature=features.reindex(x.index, fill_value=0.0),
            )
        elif algorithm == "modularity":
            return prefix + busmap_by_greedy_modularity(
                n, n_clusters_c[x.name], buses_i=x.index
            )
        else:
            raise ValueError(
                f"`algorithm` must be one of 'kmeans' or 'hac' or 'modularity'. Is {algorithm}."
            )

    compat_kws = dict(include_groups=False) if PD_GE_2_2 else {}

    return (
        n.buses.groupby(["country", "sub_network"], group_keys=False)
        .apply(busmap_for_country, **compat_kws)
        .squeeze()
        .rename("busmap")
    )


def clustering_for_n_clusters(
    n: pypsa.Network,
    busmap: pd.Series,
    aggregation_strategies: dict | None = None,
) -> pypsa.clustering.spatial.Clustering:
    if aggregation_strategies is None:
        aggregation_strategies = dict()

    line_strategies = aggregation_strategies.get("lines", dict())

    bus_strategies = aggregation_strategies.get("buses", dict())
    bus_strategies.setdefault("substation_lv", lambda x: bool(x.sum()))
    bus_strategies.setdefault("substation_off", lambda x: bool(x.sum()))

    clustering = get_clustering_from_busmap(
        n,
        busmap,
        bus_strategies=bus_strategies,
        line_strategies=line_strategies,
        custom_line_groupers=["build_year"],
    )

    return clustering


def cluster_regions(
    busmaps: tuple | list, regions: gpd.GeoDataFrame, with_country: bool = False
) -> gpd.GeoDataFrame:
    """
    Cluster regions based on busmaps and save the results to a file and to the
    network.

    Parameters
    ----------
    - busmaps (list): A list of busmaps used for clustering.
    - regions (gpd.GeoDataFrame): The regions to cluster.
    - with_country (bool): Whether to keep country column.

    Returns
    -------
    None
    """
    busmap = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])
    columns = ["name", "country", "geometry"] if with_country else ["name", "geometry"]
    regions = regions.reindex(columns=columns).set_index("name")
    regions_c = regions.dissolve(busmap)
    regions_c.index.name = "name"
    return regions_c.reset_index()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("cluster_network", clusters=60)
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params
    solver_name = snakemake.config["solving"]["solver"]["name"]

    n = pypsa.Network(snakemake.input.network)
    buses_prev, lines_prev, links_prev = len(n.buses), len(n.lines), len(n.links)

    load = (
        xr.open_dataarray(snakemake.input.load)
        .mean(dim="time")
        .to_pandas()
        .reindex(n.buses.index, fill_value=0.0)
    )

    if snakemake.wildcards.clusters == "all":
        n_clusters = len(n.buses)
    else:
        n_clusters = int(snakemake.wildcards.clusters)

    if n_clusters == len(n.buses):
        # Fast-path if no clustering is necessary
        busmap = n.buses.index.to_series()
        linemap = n.lines.index.to_series()
        clustering = pypsa.clustering.spatial.Clustering(n, busmap, linemap)
    else:
        Nyears = n.snapshot_weightings.objective.sum() / 8760

        custom_busmap = params.custom_busmap
        if custom_busmap:
            custom_busmap = pd.read_csv(
                snakemake.input.custom_busmap, index_col=0
            ).squeeze()
            custom_busmap.index = custom_busmap.index.astype(str)
            logger.info(f"Imported custom busmap from {snakemake.input.custom_busmap}")
            busmap = custom_busmap
        else:
            algorithm = params.cluster_network["algorithm"]
            features = None
            if algorithm == "hac":
                features = get_feature_data_for_hac(snakemake.input.hac_features)
                fix_country_assignment_for_hac(n)

            n.determine_network_topology()

            n_clusters_c = distribute_n_clusters_to_countries(
                n,
                n_clusters,
                load,
                focus_weights=params.focus_weights,
                solver_name=solver_name,
            )

            busmap = busmap_for_n_clusters(
                n,
                n_clusters_c,
                cluster_weights=load,
                algorithm=algorithm,
                features=features,
            )

        clustering = clustering_for_n_clusters(
            n,
            busmap,
            aggregation_strategies=params.aggregation_strategies,
        )

    nc = clustering.n

    for attr in ["busmap", "linemap"]:
        getattr(clustering, attr).to_csv(snakemake.output[attr])

    # nc.shapes = n.shapes.copy()
    for which in ["regions_onshore", "regions_offshore"]:
        regions = gpd.read_file(snakemake.input[which])
        clustered_regions = cluster_regions((clustering.busmap,), regions)
        clustered_regions.to_file(snakemake.output[which])
        # append_bus_shapes(nc, clustered_regions, type=which.split("_")[1])

    nc.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    nc.export_to_netcdf(snakemake.output.network)

    logger.info(
        f"Clustered network:\n"
        f"Buses: {buses_prev} to {len(nc.buses)}\n"
        f"Lines: {lines_prev} to {len(nc.lines)}\n"
        f"Links: {links_prev} to {len(nc.links)}"
    )
