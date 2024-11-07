# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Lifts electrical transmission network to a single 380 kV voltage layer, removes
dead-ends of the network, and reduces multi-hop HVDC connections to a single
link.

Relevant Settings
-----------------

.. code:: yaml

    clustering:
      simplify_network:
      cluster_network:
      aggregation_strategies:

    links:
        p_max_pu:


.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`electricity_cf`, :ref:`renewable_cf`,
    :ref:`lines_cf`, :ref:`links_cf`

Inputs
------

- ``resources/regions_onshore.geojson``: confer :ref:`busregions`
- ``resources/regions_offshore.geojson``: confer :ref:`busregions`
- ``networks/base.nc``

Outputs
-------

- ``resources/regions_onshore_base.geojson``:

    .. image:: img/regions_onshore_base_s.png
            :scale: 33 %

- ``resources/regions_offshore_base.geojson``:

    .. image:: img/regions_offshore_base_s  .png
            :scale: 33 %

- ``resources/busmap_base_s.csv``: Mapping of buses from ``networks/base.nc`` to ``networks/base_s.nc``;
- ``networks/base.nc``:

    .. image:: img/base_s.png
        :scale: 33 %

Description
-----------

The rule :mod:`simplify_network` does up to three things:

1. Create an equivalent transmission network in which all voltage levels are mapped to the 380 kV level by the function ``simplify_network(...)``.

2. DC only sub-networks that are connected at only two buses to the AC network are reduced to a single representative link in the function ``simplify_links(...)``.

3. Stub lines and links, i.e. dead-ends of the network, are sequentially removed from the network in the function ``remove_stubs(...)``.
"""

import logging
from functools import reduce
from typing import Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import scipy as sp
from _helpers import configure_logging, set_scenario_config
from base_network import append_bus_shapes
from cluster_network import cluster_regions
from pypsa.clustering.spatial import busmap_by_stubs, get_clustering_from_busmap
from scipy.sparse.csgraph import connected_components, dijkstra

logger = logging.getLogger(__name__)


def simplify_network_to_380(
    n: pypsa.Network, linetype_380: str
) -> Tuple[pypsa.Network, pd.Series]:
    """
    Fix all lines to a voltage level of 380 kV and remove all transformers.

    The function preserves the transmission capacity for each line while
    updating its voltage level, line type and number of parallel bundles
    (num_parallel).

    Transformers are removed and connected components are moved from
    their starting bus to their ending bus. The corresponding starting
    buses are removed as well.
    """
    logger.info("Mapping all network lines onto a single 380kV layer")

    n.buses["v_nom"] = 380.0

    n.lines["type"] = linetype_380
    n.lines["v_nom"] = 380
    n.lines["i_nom"] = n.line_types.i_nom[linetype_380]
    n.lines["num_parallel"] = n.lines.eval("s_nom / (sqrt(3) * v_nom * i_nom)")

    trafo_map = pd.Series(n.transformers.bus1.values, n.transformers.bus0.values)
    trafo_map = trafo_map[~trafo_map.index.duplicated(keep="first")]
    while (several_trafo_b := trafo_map.isin(trafo_map.index)).any():
        trafo_map[several_trafo_b] = trafo_map[several_trafo_b].map(trafo_map)
    missing_buses_i = n.buses.index.difference(trafo_map.index)
    missing = pd.Series(missing_buses_i, missing_buses_i)
    trafo_map = pd.concat([trafo_map, missing])

    for c in n.one_port_components | n.branch_components:
        df = n.df(c)
        for col in df.columns:
            if col.startswith("bus"):
                df[col] = df[col].map(trafo_map)

    n.remove("Transformer", n.transformers.index)
    n.remove("Bus", n.buses.index.difference(trafo_map))

    return n, trafo_map


def _remove_clustered_buses_and_branches(n: pypsa.Network, busmap: pd.Series) -> None:
    buses_to_del = n.buses.index.difference(busmap)
    n.remove("Bus", buses_to_del)
    for c in n.branch_components:
        df = n.df(c)
        n.remove(c, df.index[df.bus0.isin(buses_to_del) | df.bus1.isin(buses_to_del)])


def simplify_links(
    n: pypsa.Network, p_max_pu: int | float
) -> Tuple[pypsa.Network, pd.Series]:
    ## Complex multi-node links are folded into end-points
    logger.info("Simplifying connected link components")

    if n.links.empty:
        return n, n.buses.index.to_series()

    # Determine connected link components, ignore all links but DC
    adjacency_matrix = n.adjacency_matrix(
        branch_components=["Link"],
        weights=dict(Link=(n.links.carrier == "DC").astype(float)),
    )

    _, labels = connected_components(adjacency_matrix, directed=False)
    labels = pd.Series(labels, n.buses.index)

    # Only span graph over the DC link components
    G = n.graph(branch_components=["Link"])

    def split_links(nodes, added_supernodes):
        nodes = frozenset(nodes)

        seen = set()

        # Supernodes are endpoints of links, identified by having lass then two neighbours or being an AC Bus
        # An example for the latter is if two different links are connected to the same AC bus.
        supernodes = {
            m
            for m in nodes
            if (
                (len(G.adj[m]) < 2 or (set(G.adj[m]) - nodes))
                or (n.buses.loc[m, "carrier"] == "AC")
                or (m in added_supernodes)
            )
        }

        for u in supernodes:
            for m, ls in G.adj[u].items():
                if m not in nodes or m in seen:
                    continue

                buses = [u, m]
                links = [list(ls)]  # [name for name in ls]]

                while m not in (supernodes | seen):
                    seen.add(m)
                    for m2, ls in G.adj[m].items():
                        if m2 in seen or m2 == u:
                            continue
                        buses.append(m2)
                        links.append(list(ls))  # [name for name in ls])
                        break
                    else:
                        # stub
                        break
                    m = m2
                if m != u:
                    yield pd.Index((u, m)), buses, links
            seen.add(u)

    busmap = n.buses.index.to_series()

    node_corsica = find_closest_bus(
        n,
        x=9.44802,
        y=42.52842,
        tol=2000,  # Tolerance needed to only return the bus if the region is actually modelled
    )

    added_supernodes = []
    if node_corsica is not None:
        added_supernodes.append(node_corsica)

    for lbl in labels.value_counts().loc[lambda s: s > 2].index:
        for b, buses, links in split_links(
            labels.index[labels == lbl], added_supernodes
        ):
            if len(buses) <= 2:
                continue

            logger.debug("nodes = {}".format(labels.index[labels == lbl]))
            logger.debug("b = {}\nbuses = {}\nlinks = {}".format(b, buses, links))

            m = sp.spatial.distance_matrix(
                n.buses.loc[b, ["x", "y"]], n.buses.loc[buses[1:-1], ["x", "y"]]
            )
            busmap.loc[buses] = b[np.r_[0, m.argmin(axis=0), 1]]

            all_links = [i for _, i in sum(links, [])]

            lengths = n.links.loc[all_links, "length"]
            name = lengths.idxmax() + "+{}".format(len(links) - 1)
            params = dict(
                carrier="DC",
                bus0=b[0],
                bus1=b[1],
                length=sum(
                    n.links.loc[[i for _, i in l], "length"].mean() for l in links
                ),
                p_nom=min(n.links.loc[[i for _, i in l], "p_nom"].sum() for l in links),
                underwater_fraction=sum(
                    lengths
                    / lengths.sum()
                    * n.links.loc[all_links, "underwater_fraction"]
                ),
                p_max_pu=p_max_pu,
                p_min_pu=-p_max_pu,
                underground=False,
                under_construction=False,
            )

            logger.info(
                "Joining the links {} connecting the buses {} to simple link {}".format(
                    ", ".join(all_links), ", ".join(buses), name
                )
            )

            n.remove("Link", all_links)

            static_attrs = n.components["Link"]["attrs"].loc[lambda df: df.static]
            for attr, default in static_attrs.default.items():
                params.setdefault(attr, default)
            n.links.loc[name] = pd.Series(params)

            # n.add("Link", name, **params)

    logger.debug("Collecting all components using the busmap")

    _remove_clustered_buses_and_branches(n, busmap)

    # Change carrier type of all added super_nodes to "AC"
    n.buses.loc[added_supernodes, "carrier"] = "AC"

    return n, busmap


def remove_stubs(
    n: pypsa.Network, simplify_network: dict
) -> Tuple[pypsa.Network, pd.Series]:
    logger.info("Removing stubs")

    across_borders = simplify_network["remove_stubs_across_borders"]
    matching_attrs = [] if across_borders else ["country"]
    busmap = busmap_by_stubs(n, matching_attrs)

    _remove_clustered_buses_and_branches(n, busmap)

    return n, busmap


def aggregate_to_substations(
    n: pypsa.Network,
    buses_i: pd.Index | list,
    aggregation_strategies: dict | None = None,
) -> Tuple[pypsa.Network, pd.Series]:
    # can be used to aggregate a selection of buses to electrically closest neighbors
    logger.info("Aggregating buses to substations")
    if aggregation_strategies is None:
        aggregation_strategies = dict()

    weight = pd.concat(
        {
            "Line": n.lines.length / n.lines.s_nom.clip(1e-3),
            "Link": n.links.length / n.links.p_nom.clip(1e-3),
        }
    )

    adj = n.adjacency_matrix(branch_components=["Line", "Link"], weights=weight)

    bus_indexer = n.buses.index.get_indexer(buses_i)
    dist = pd.DataFrame(
        dijkstra(adj, directed=False, indices=bus_indexer), buses_i, n.buses.index
    )

    dist[buses_i] = (
        np.inf
    )  # bus in buses_i should not be assigned to different bus in buses_i

    for c in n.buses.country.unique():
        incountry_b = n.buses.country == c
        dist.loc[incountry_b, ~incountry_b] = np.inf

    busmap = n.buses.index.to_series()
    busmap.loc[buses_i] = dist.idxmin(1)

    line_strategies = aggregation_strategies.get("lines", dict())

    bus_strategies = aggregation_strategies.get("buses", dict())
    bus_strategies.setdefault("substation_lv", lambda x: bool(x.sum()))
    bus_strategies.setdefault("substation_off", lambda x: bool(x.sum()))

    clustering = get_clustering_from_busmap(
        n,
        busmap,
        line_length_factor=1.0,
        bus_strategies=bus_strategies,
        line_strategies=line_strategies,
    )
    return clustering.n, busmap


def find_closest_bus(n, x, y, tol=2000):
    """
    Find the index of the closest bus to the given coordinates within a specified tolerance.
    Parameters:
        n (pypsa.Network): The network object.
        x (float): The x-coordinate (longitude) of the target location.
        y (float): The y-coordinate (latitude) of the target location.
        tol (float): The distance tolerance in meters. Default is 2000 meters.

    Returns:
        int: The index of the closest bus to the target location within the tolerance.
             Returns None if no bus is within the tolerance.
    """
    # Conversion factors
    meters_per_degree_lat = 111139  # Meters per degree of latitude
    meters_per_degree_lon = 111139 * np.cos(
        np.radians(y)
    )  # Meters per degree of longitude at the given latitude

    x0 = np.array(n.buses.x)
    y0 = np.array(n.buses.y)

    # Calculate distances in meters
    dist = np.sqrt(
        ((x - x0) * meters_per_degree_lon) ** 2
        + ((y - y0) * meters_per_degree_lat) ** 2
    )

    # Find the closest bus within the tolerance
    min_dist = dist.min()
    if min_dist <= tol:
        return n.buses.index[dist.argmin()]
    else:
        return None


def remove_converters(n: pypsa.Network) -> pypsa.Network:
    """
    Remove all converters from the network and remap all buses that were originally connected to the
    converter to the connected AC bus. Preparation step before simplifying links.

    Parameters:
        n (pypsa.Network): The network object.

    Returns:
        n (pypsa.Network): The network object with all converters removed.
    """
    # Extract converters
    converters = n.links.query("carrier == ''")[["bus0", "bus1"]]
    converters["bus0_carrier"] = converters["bus0"].map(n.buses.carrier)
    converters["bus1_carrier"] = converters["bus1"].map(n.buses.carrier)

    converters["ac_bus"] = converters.apply(
        lambda x: x["bus1"] if x["bus1_carrier"] == "AC" else x["bus0"], axis=1
    )

    converters["dc_bus"] = converters.apply(
        lambda x: x["bus1"] if x["bus1_carrier"] == "DC" else x["bus0"], axis=1
    )

    # Dictionary for remapping
    dict_dc_to_ac = dict(zip(converters["dc_bus"], converters["ac_bus"]))

    # Remap all buses that were originally connected to the converter to the connected AC bus
    n.links["bus0"] = n.links["bus0"].replace(dict_dc_to_ac)
    n.links["bus1"] = n.links["bus1"].replace(dict_dc_to_ac)

    # Remove all converters from network.links and associated dc buses from network.buses
    n.links = n.links.loc[~n.links.index.isin(converters.index)]
    n.buses = n.buses.loc[~n.buses.index.isin(converters["dc_bus"])]

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("simplify_network")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params

    n = pypsa.Network(snakemake.input.network)
    Nyears = n.snapshot_weightings.objective.sum() / 8760
    buses_prev, lines_prev, links_prev = len(n.buses), len(n.lines), len(n.links)

    linetype_380 = snakemake.config["lines"]["types"][380]
    n, trafo_map = simplify_network_to_380(n, linetype_380)

    n = remove_converters(n)

    n, simplify_links_map = simplify_links(n, params.p_max_pu)

    busmaps = [trafo_map, simplify_links_map]

    if params.simplify_network["remove_stubs"]:
        n, stub_map = remove_stubs(n, params.simplify_network)
        busmaps.append(stub_map)

    substations_i = n.buses.query("substation_lv or substation_off").index

    # some entries in n.buses are not updated in previous functions, therefore can be wrong. as they are not needed
    # and are lost when clustering (for example with the simpl wildcard), we remove them for consistency:
    remove = [
        "symbol",
        "tags",
        "under_construction",
        "onshore_bus",
        "geometry",
        "underground",
        "project_status",
    ]
    n.buses.drop(remove, axis=1, inplace=True, errors="ignore")
    n.lines.drop(remove, axis=1, errors="ignore", inplace=True)

    if params.simplify_network["to_substations"]:
        n, substation_map = aggregate_to_substations(
            n, substations_i, params.aggregation_strategies
        )
        busmaps.append(substation_map)

    # all buses without shapes need to be clustered to their closest neighbor for HAC
    if params.cluster_network["algorithm"] == "hac":
        buses_i = list(n.buses.index.difference(n.shapes.idx))
        logger.info(
            "Preparing for HAC-Clustering. "
            f"Aggregating {len(buses_i)} buses without Voronoi shapes to closest neighbor."
        )
        n, busmap_hac = aggregate_to_substations(
            n, buses_i, params.aggregation_strategies
        )
        busmaps.append(busmap_hac)

    busmap_s = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])
    busmap_s.to_csv(snakemake.output.busmap)

    for which in ["regions_onshore", "regions_offshore"]:
        regions = gpd.read_file(snakemake.input[which])
        clustered_regions = cluster_regions(busmaps, regions, with_country=True)
        clustered_regions.to_file(snakemake.output[which])
        # append_bus_shapes(n, clustered_regions, type=which.split("_")[1])

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)

    logger.info(
        f"Simplified network:\n"
        f"Buses: {buses_prev} to {len(n.buses)}\n"
        f"Lines: {lines_prev} to {len(n.lines)}\n"
        f"Links: {links_prev} to {len(n.links)}"
    )
