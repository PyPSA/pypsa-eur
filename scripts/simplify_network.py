# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

# coding: utf-8
"""
Lifts electrical transmission network to a single 380 kV voltage layer,
removes dead-ends of the network,
and reduces multi-hop HVDC connections to a single link.

Relevant Settings
-----------------

.. code:: yaml

    costs:
        USD2013_to_EUR2013:
        discountrate:
        marginal_cost:
        capital_cost:

    electricity:
        max_hours:

    renewables: (keys)
        {technology}:
            potential:

    lines:
        length_factor:

    links:
        p_max_pu:

    solving:
        solver:
            name:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`costs_cf`, :ref:`electricity_cf`, :ref:`renewable_cf`,
    :ref:`lines_cf`, :ref:`links_cf`, :ref:`solving_cf`

Inputs
------

- ``data/costs.csv``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``resources/regions_onshore.geojson``: confer :ref:`busregions`
- ``resources/regions_offshore.geojson``: confer :ref:`busregions`
- ``networks/elec.nc``: confer :ref:`electricity`

Outputs
-------

- ``resources/regions_onshore_elec_s{simpl}.geojson``:

    .. image:: ../img/regions_onshore_elec_s.png
            :scale: 33 %

- ``resources/regions_offshore_elec_s{simpl}.geojson``:

    .. image:: ../img/regions_offshore_elec_s  .png
            :scale: 33 %

- ``resources/busmap_elec_s{simpl}.csv``: Mapping of buses from ``networks/elec.nc`` to ``networks/elec_s{simpl}.nc``;
- ``networks/elec_s{simpl}.nc``:

    .. image:: ../img/elec_s.png
        :scale: 33 %

Description
-----------

The rule :mod:`simplify_network` does up to four things:

1. Create an equivalent transmission network in which all voltage levels are mapped to the 380 kV level by the function ``simplify_network(...)``.

2. DC only sub-networks that are connected at only two buses to the AC network are reduced to a single representative link in the function ``simplify_links(...)``. The components attached to buses in between are moved to the nearest endpoint. The grid connection cost of offshore wind generators are added to the captial costs of the generator.

3. Stub lines and links, i.e. dead-ends of the network, are sequentially removed from the network in the function ``remove_stubs(...)``. Components are moved along.

4. Optionally, if an integer were provided for the wildcard ``{simpl}`` (e.g. ``networks/elec_s500.nc``), the network is clustered to this number of clusters with the routines from the ``cluster_network`` rule with the function ``cluster_network.cluster(...)``. This step is usually skipped!
"""

import logging
from _helpers import configure_logging

from cluster_network import clustering_for_n_clusters, cluster_regions
from add_electricity import load_costs

import pandas as pd
import numpy as np
import scipy as sp
from scipy.sparse.csgraph import connected_components, dijkstra

from functools import reduce

import pypsa
from pypsa.io import import_components_from_dataframe, import_series_from_dataframe
from pypsa.networkclustering import stubs_clustering, get_clustering_from_busmap

logger = logging.getLogger(__name__)


def simplify_network_to_380(n):
    ## All goes to v_nom == 380
    logger.info("Mapping all network lines onto a single 380kV layer")

    n.buses['v_nom'] = 380.

    linetype_380, = n.lines.loc[n.lines.v_nom == 380., 'type'].unique()
    lines_v_nom_b = n.lines.v_nom != 380.
    n.lines.loc[lines_v_nom_b, 'num_parallel'] *= (n.lines.loc[lines_v_nom_b, 'v_nom'] / 380.)**2
    n.lines.loc[lines_v_nom_b, 'v_nom'] = 380.
    n.lines.loc[lines_v_nom_b, 'type'] = linetype_380
    n.lines.loc[lines_v_nom_b, 's_nom'] = (
        np.sqrt(3) * n.lines['type'].map(n.line_types.i_nom) *
        n.lines.bus0.map(n.buses.v_nom) * n.lines.num_parallel
    )

    # Replace transformers by lines
    trafo_map = pd.Series(n.transformers.bus1.values, index=n.transformers.bus0.values)
    trafo_map = trafo_map[~trafo_map.index.duplicated(keep='first')]
    several_trafo_b = trafo_map.isin(trafo_map.index)
    trafo_map.loc[several_trafo_b] = trafo_map.loc[several_trafo_b].map(trafo_map)
    missing_buses_i = n.buses.index.difference(trafo_map.index)
    trafo_map = trafo_map.append(pd.Series(missing_buses_i, missing_buses_i))

    for c in n.one_port_components|n.branch_components:
        df = n.df(c)
        for col in df.columns:
            if col.startswith('bus'):
                df[col] = df[col].map(trafo_map)

    n.mremove("Transformer", n.transformers.index)
    n.mremove("Bus", n.buses.index.difference(trafo_map))

    return n, trafo_map

def simplify_links(n):
    logger.info("Simplifying connected link components")

    if n.links.empty:
        return n, n.buses.index.to_series()

    busmap = n.buses.index.to_series()

    n.buses.carrier = 'AC'
    
    c1_buses = []
    for m in n.buses.index:
        links = n.links[(n.links.bus0==m) | (n.links.bus1==m)].index
        lines = n.lines[(n.lines.bus0==m) | (n.lines.bus1==m)].index
        if (len(links) > 1) & (len(lines) == 0): # if == 1 -> remove_stubs
            c1_buses.append(m)
            link = n.links.loc[links].index[n.links.loc[links].length.argmin()]
            candidates = set(n.links.loc[link][['bus0','bus1']])-set(c1_buses)
            if len(candidates) == 0:
                link = n.links.loc[links].index[n.links.loc[links].length.argmax()]
                candidates = set(n.links.loc[link][['bus0','bus1']])-set(c1_buses)
            print(m, c1_buses, len(links))
            bus = list(candidates)[0]
            if n.buses.loc[m].country == n.buses.loc[bus].country:
                busmap[m] = bus

    busmap = busmap.map(busmap)
    clustering = get_clustering_from_busmap(n, busmap, bus_strategies=dict(country=lambda x: x.value_counts().index[0]))

    return clustering.network, busmap
                
            

def remove_stubs(n):
    logger.info("Removing stubs")

    n.buses.carrier = 'AC'
    clustering = stubs_clustering(n)
    
    return clustering.network, clustering.busmap


def cluster(n, n_clusters):
    logger.info(f"Clustering to {n_clusters} buses")

    focus_weights = snakemake.config.get('focus_weights', None)
    
    renewable_carriers = pd.Index([tech
                                    for tech in n.generators.carrier.unique()
                                    if tech.split('-', 2)[0] in snakemake.config['renewable']])
    def consense(x):
        v = x.iat[0]
        assert ((x == v).all() or x.isnull().all()), (
            "The `potential` configuration option must agree for all renewable carriers, for now!"
        )
        return v
    potential_mode = (consense(pd.Series([snakemake.config['renewable'][tech]['potential']
                                            for tech in renewable_carriers]))
                        if len(renewable_carriers) > 0 else 'conservative')
    clustering = clustering_for_n_clusters(n, n_clusters, custom_busmap=False, potential_mode=potential_mode,
                                           solver_name=snakemake.config['solving']['solver']['name'],
                                           focus_weights=focus_weights)

    return clustering.network, clustering.busmap


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('simplify_network', simpl='', network='elec')
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    n, trafo_map = simplify_network_to_380(n)

    n, simplify_links_map = simplify_links(n)

    n, stub_map = remove_stubs(n)

    busmaps = [trafo_map, simplify_links_map, stub_map]

    if snakemake.wildcards.simpl:
        n, cluster_map = cluster(n, int(snakemake.wildcards.simpl))
        busmaps.append(cluster_map)

    n.export_to_netcdf(snakemake.output.network)

    busmap_s = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])
    busmap_s.to_csv(snakemake.output.busmap)

    cluster_regions(busmaps, snakemake.input, snakemake.output)
