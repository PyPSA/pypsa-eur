# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Creates the network topology from a `ENTSO-E map extract.
<https://github.com/PyPSA/GridKit/tree/master/entsoe>`_ (March 2022)
or `OpenStreetMap data <https://www.openstreetmap.org/>`_ (Aug 2024)
as a PyPSA
network.

Relevant Settings
-----------------

.. code:: yaml

    countries:

    electricity:
        voltages:

    lines:
        types:
        s_max_pu:
        under_construction:

    links:
        p_max_pu:
        under_construction:

    transformers:
        x:
        s_nom:
        type:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`snapshots_cf`, :ref:`toplevel_cf`, :ref:`electricity_cf`, :ref:`load_cf`,
    :ref:`lines_cf`, :ref:`links_cf`, :ref:`transformers_cf`

Inputs
------

- ``data/entsoegridkit``:  Extract from the geographical vector data of the online `ENTSO-E Interactive Map <https://www.entsoe.eu/data/map/>`_ by the `GridKit <https://github.com/martacki/gridkit>`_ toolkit dating back to March 2022.
- ``data/parameter_corrections.yaml``: Corrections for ``data/entsoegridkit``
- ``data/links_p_nom.csv``: confer :ref:`links`
- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``resources/offshore_shapes.geojson``: confer :ref:`shapes`
- ``resources/europe_shape.geojson``: confer :ref:`shapes`

Outputs
-------

- ``networks/base.nc``

    .. image:: img/base.png
        :scale: 33 %

- ``resources/regions_onshore.geojson``:

    .. image:: img/regions_onshore.png
        :scale: 33 %

- ``resources/regions_offshore.geojson``:

    .. image:: img/regions_offshore.png
        :scale: 33 %

Description
-----------
Creates the network topology from an ENTSO-E map extract, and create Voronoi shapes for each bus representing both onshore and offshore regions.
"""

import logging
import warnings
from itertools import product

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pypsa
import shapely
import shapely.prepared
import shapely.wkt
import yaml
from _helpers import REGION_COLS, configure_logging, get_snapshots, set_scenario_config
from packaging.version import Version, parse
from scipy.sparse import csgraph
from scipy.spatial import KDTree
from shapely.geometry import LineString, Point

PD_GE_2_2 = parse(pd.__version__) >= Version("2.2")

logger = logging.getLogger(__name__)


def _get_oid(df):
    if "tags" in df.columns:
        return df.tags.str.extract(r'"oid"=>"(\d+)"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def _get_country(df):
    if "tags" in df.columns:
        return df.tags.str.extract('"country"=>"([A-Z]{2})"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def _find_closest_links(links, new_links, distance_upper_bound=1.5):
    treecoords = np.asarray(
        [
            np.asarray(shapely.wkt.loads(s).coords)[[0, -1]].flatten()
            for s in links.geometry
        ]
    )
    querycoords = np.vstack(
        [new_links[["x1", "y1", "x2", "y2"]], new_links[["x2", "y2", "x1", "y1"]]]
    )
    tree = KDTree(treecoords)
    dist, ind = tree.query(querycoords, distance_upper_bound=distance_upper_bound)
    found_b = ind < len(links)
    found_i = np.arange(len(new_links) * 2)[found_b] % len(new_links)
    return (
        pd.DataFrame(
            dict(D=dist[found_b], i=links.index[ind[found_b] % len(links)]),
            index=new_links.index[found_i],
        )
        .sort_values(by="D")[lambda ds: ~ds.index.duplicated(keep="first")]
        .sort_index()["i"]
    )


def _load_buses(buses, europe_shape, countries, config):
    buses = (
        pd.read_csv(
            buses,
            quotechar="'",
            true_values=["t"],
            false_values=["f"],
            dtype=dict(bus_id="str"),
        )
        .set_index("bus_id")
        .rename(columns=dict(voltage="v_nom"))
    )

    if "station_id" in buses.columns:
        buses.drop("station_id", axis=1, inplace=True)

    buses["carrier"] = buses.pop("dc").map({True: "DC", False: "AC"})
    buses["under_construction"] = buses.under_construction.where(
        lambda s: s.notnull(), False
    ).astype(bool)
    europe_shape = gpd.read_file(europe_shape).loc[0, "geometry"]
    europe_shape_prepped = shapely.prepared.prep(europe_shape)
    buses_in_europe_b = buses[["x", "y"]].apply(
        lambda p: europe_shape_prepped.contains(Point(p)), axis=1
    )

    buses_in_countries_b = (
        buses.country.isin(countries)
        if "country" in buses
        else pd.Series(True, buses.index)
    )

    v_nom_min = min(config["electricity"]["voltages"])
    v_nom_max = max(config["electricity"]["voltages"])

    buses_with_v_nom_to_keep_b = (
        (v_nom_min <= buses.v_nom) & (buses.v_nom <= v_nom_max)
        | (buses.v_nom.isnull())
        | (
            buses.carrier == "DC"
        )  # Keeping all DC buses from the input dataset independent of voltage (e.g. 150 kV connections)
    )

    logger.info(f"Removing buses outside of range AC {v_nom_min} - {v_nom_max} V")
    return pd.DataFrame(
        buses.loc[buses_in_europe_b & buses_in_countries_b & buses_with_v_nom_to_keep_b]
    )


def _load_transformers(buses, transformers):
    transformers = pd.read_csv(
        transformers,
        quotechar="'",
        true_values=["t"],
        false_values=["f"],
        dtype=dict(transformer_id="str", bus0="str", bus1="str"),
    ).set_index("transformer_id")

    transformers = _remove_dangling_branches(transformers, buses)

    return transformers


def _load_converters_from_eg(buses, converters):
    converters = pd.read_csv(
        converters,
        quotechar="'",
        true_values=["t"],
        false_values=["f"],
        dtype=dict(converter_id="str", bus0="str", bus1="str"),
    ).set_index("converter_id")

    converters = _remove_dangling_branches(converters, buses)

    converters["carrier"] = "B2B"

    return converters


def _load_converters_from_osm(buses, converters):
    converters = pd.read_csv(
        converters,
        quotechar="'",
        true_values=["t"],
        false_values=["f"],
        dtype=dict(converter_id="str", bus0="str", bus1="str"),
    ).set_index("converter_id")

    converters = _remove_dangling_branches(converters, buses)

    converters["carrier"] = ""

    return converters


def _load_links_from_eg(buses, links):
    links = pd.read_csv(
        links,
        quotechar="'",
        true_values=["t"],
        false_values=["f"],
        dtype=dict(link_id="str", bus0="str", bus1="str", under_construction="bool"),
    ).set_index("link_id")

    links["length"] /= 1e3

    # Skagerrak Link is connected to 132kV bus which is removed in _load_buses.
    # Connect to neighboring 380kV bus
    links.loc[links.bus1 == "6396", "bus1"] = "6398"

    links = _remove_dangling_branches(links, buses)

    # Add DC line parameters
    links["carrier"] = "DC"

    return links


def _load_links_from_osm(buses, links):
    links = pd.read_csv(
        links,
        quotechar="'",
        true_values=["t"],
        false_values=["f"],
        dtype=dict(
            link_id="str",
            bus0="str",
            bus1="str",
            voltage="int",
            p_nom="float",
        ),
    ).set_index("link_id")

    links["length"] /= 1e3

    links = _remove_dangling_branches(links, buses)

    # Add DC line parameters
    links["carrier"] = "DC"

    return links


def _load_lines(buses, lines):
    lines = (
        pd.read_csv(
            lines,
            quotechar="'",
            true_values=["t"],
            false_values=["f"],
            dtype=dict(
                line_id="str",
                bus0="str",
                bus1="str",
                underground="bool",
                under_construction="bool",
            ),
        )
        .set_index("line_id")
        .rename(columns=dict(voltage="v_nom", circuits="num_parallel"))
    )

    lines["length"] /= 1e3

    lines["carrier"] = "AC"
    lines = _remove_dangling_branches(lines, buses)

    return lines


def _apply_parameter_corrections(n, parameter_corrections):
    with open(parameter_corrections) as f:
        corrections = yaml.safe_load(f)

    if corrections is None:
        return

    for component, attrs in corrections.items():
        df = n.df(component)
        oid = _get_oid(df)
        if attrs is None:
            continue

        for attr, repls in attrs.items():
            for i, r in repls.items():
                if i == "oid":
                    r = oid.map(repls["oid"]).dropna()
                elif i == "index":
                    r = pd.Series(repls["index"])
                else:
                    raise NotImplementedError()
                inds = r.index.intersection(df.index)
                df.loc[inds, attr] = r[inds].astype(df[attr].dtype)


def _reconnect_crimea(lines):
    logger.info("Reconnecting Crimea to the Ukrainian grid.")
    lines_to_crimea = pd.DataFrame(
        {
            "bus0": ["3065", "3181", "3181"],
            "bus1": ["3057", "3055", "3057"],
            "v_nom": [300, 300, 300],
            "num_parallel": [1, 1, 1],
            "length": [140, 120, 140],
            "carrier": ["AC", "AC", "AC"],
            "underground": [False, False, False],
            "under_construction": [False, False, False],
        },
        index=["Melitopol", "Liubymivka left", "Luibymivka right"],
    )

    return pd.concat([lines, lines_to_crimea])


def _set_electrical_parameters_lines_eg(lines, config):
    v_noms = config["electricity"]["voltages"]
    linetypes = config["lines"]["types"]

    for v_nom in v_noms:
        lines.loc[lines["v_nom"] == v_nom, "type"] = linetypes[v_nom]

    lines["s_max_pu"] = config["lines"]["s_max_pu"]

    return lines


def _set_electrical_parameters_lines_osm(lines, config):
    if lines.empty:
        lines["type"] = []
        return lines

    v_noms = config["electricity"]["voltages"]
    linetypes = _get_linetypes_config(config["lines"]["types"], v_noms)

    lines["carrier"] = "AC"
    lines["dc"] = False

    lines.loc[:, "type"] = lines.v_nom.apply(
        lambda x: _get_linetype_by_voltage(x, linetypes)
    )

    lines["s_max_pu"] = config["lines"]["s_max_pu"]

    return lines


def _set_lines_s_nom_from_linetypes(n):
    n.lines["s_nom"] = (
        np.sqrt(3)
        * n.lines["type"].map(n.line_types.i_nom)
        * n.lines["v_nom"]
        * n.lines["num_parallel"]
    )


def _set_electrical_parameters_links_eg(links, config, links_p_nom):
    if links.empty:
        return links

    p_max_pu = config["links"].get("p_max_pu", 1.0)
    links["p_max_pu"] = p_max_pu
    links["p_min_pu"] = -p_max_pu

    links_p_nom = pd.read_csv(links_p_nom)

    # filter links that are not in operation anymore
    removed_b = links_p_nom.Remarks.str.contains("Shut down|Replaced", na=False)
    links_p_nom = links_p_nom[~removed_b]

    # find closest link for all links in links_p_nom
    links_p_nom["j"] = _find_closest_links(links, links_p_nom)

    links_p_nom = links_p_nom.groupby(["j"], as_index=False).agg({"Power (MW)": "sum"})

    p_nom = links_p_nom.dropna(subset=["j"]).set_index("j")["Power (MW)"]

    # Don't update p_nom if it's already set
    p_nom_unset = (
        p_nom.drop(links.index[links.p_nom.notnull()], errors="ignore")
        if "p_nom" in links
        else p_nom
    )
    links.loc[p_nom_unset.index, "p_nom"] = p_nom_unset

    return links


def _set_electrical_parameters_links_osm(links, config):
    if links.empty:
        return links

    p_max_pu = config["links"].get("p_max_pu", 1.0)
    links["p_max_pu"] = p_max_pu
    links["p_min_pu"] = -p_max_pu
    links["carrier"] = "DC"
    links["dc"] = True

    return links


def _set_electrical_parameters_converters(converters, config):
    p_max_pu = config["links"].get("p_max_pu", 1.0)
    converters["p_max_pu"] = p_max_pu
    converters["p_min_pu"] = -p_max_pu

    # if column "p_nom" does not exist, set to 2000
    if "p_nom" not in converters:
        converters["p_nom"] = 2000

    # Converters are combined with links
    converters["under_construction"] = False
    converters["underground"] = False

    return converters


def _set_electrical_parameters_transformers(transformers, config):
    config = config["transformers"]

    ## Add transformer parameters
    transformers["x"] = config.get("x", 0.1)
    if "s_nom" not in transformers:
        transformers["s_nom"] = config.get("s_nom", 2000)
    transformers["type"] = config.get("type", "")

    return transformers


def _remove_dangling_branches(branches, buses):
    return pd.DataFrame(
        branches.loc[branches.bus0.isin(buses.index) & branches.bus1.isin(buses.index)]
    )


def _remove_unconnected_components(network, threshold=6):
    _, labels = csgraph.connected_components(network.adjacency_matrix(), directed=False)
    component = pd.Series(labels, index=network.buses.index)

    component_sizes = component.value_counts()
    components_to_remove = component_sizes.loc[component_sizes < threshold]

    logger.info(
        f"Removing {len(components_to_remove)} unconnected network components with less than {components_to_remove.max()} buses. In total {components_to_remove.sum()} buses."
    )

    return network[component == component_sizes.index[0]]


def _set_countries_and_substations(n, config, country_shapes, offshore_shapes):
    buses = n.buses

    def buses_in_shape(shape):
        shape = shapely.prepared.prep(shape)
        return pd.Series(
            np.fromiter(
                (
                    shape.contains(Point(x, y))
                    for x, y in buses.loc[:, ["x", "y"]].values
                ),
                dtype=bool,
                count=len(buses),
            ),
            index=buses.index,
        )

    countries = config["countries"]
    country_shapes = gpd.read_file(country_shapes).set_index("name")["geometry"]
    # reindexing necessary for supporting empty geo-dataframes
    offshore_shapes = gpd.read_file(offshore_shapes)
    offshore_shapes = offshore_shapes.reindex(columns=["name", "geometry"]).set_index(
        "name"
    )["geometry"]
    substation_b = buses["symbol"].str.contains(
        "substation|converter station", case=False
    )

    def prefer_voltage(x, which):
        index = x.index
        if len(index) == 1:
            return pd.Series(index, index)
        key = (
            x.index[0]
            if x["v_nom"].isnull().all()
            else getattr(x["v_nom"], "idx" + which)()
        )
        return pd.Series(key, index)

    compat_kws = dict(include_groups=False) if PD_GE_2_2 else {}
    gb = buses.loc[substation_b].groupby(
        ["x", "y"], as_index=False, group_keys=False, sort=False
    )
    bus_map_low = gb.apply(prefer_voltage, "min", **compat_kws)
    lv_b = (bus_map_low == bus_map_low.index).reindex(buses.index, fill_value=False)
    bus_map_high = gb.apply(prefer_voltage, "max", **compat_kws)
    hv_b = (bus_map_high == bus_map_high.index).reindex(buses.index, fill_value=False)

    onshore_b = pd.Series(False, buses.index)
    offshore_b = pd.Series(False, buses.index)

    for country in countries:
        onshore_shape = country_shapes[country]
        onshore_country_b = buses_in_shape(onshore_shape)
        onshore_b |= onshore_country_b

        buses.loc[onshore_country_b, "country"] = country

        if country not in offshore_shapes.index:
            continue
        offshore_country_b = buses_in_shape(offshore_shapes[country])
        offshore_b |= offshore_country_b

        buses.loc[offshore_country_b, "country"] = country

    # Only accept buses as low-voltage substations (where load is attached), if
    # they have at least one connection which is not under_construction
    has_connections_b = pd.Series(False, index=buses.index)
    for b, df in product(("bus0", "bus1"), (n.lines, n.links)):
        has_connections_b |= ~df.groupby(b).under_construction.min()

    buses["onshore_bus"] = onshore_b
    buses["substation_lv"] = (
        lv_b & onshore_b & (~buses["under_construction"]) & has_connections_b
    )
    buses["substation_off"] = (offshore_b | (hv_b & onshore_b)) & (
        ~buses["under_construction"]
    )

    c_nan_b = buses.country.fillna("na") == "na"
    if c_nan_b.sum() > 0:
        c_tag = _get_country(buses.loc[c_nan_b])
        c_tag.loc[~c_tag.isin(countries)] = np.nan
        n.buses.loc[c_nan_b, "country"] = c_tag

        c_tag_nan_b = n.buses.country.isnull()

        # Nearest country in path length defines country of still homeless buses
        # Work-around until commit 705119 lands in pypsa release
        n.transformers["length"] = 0.0
        graph = n.graph(weight="length")
        n.transformers.drop("length", axis=1, inplace=True)

        for b in n.buses.index[c_tag_nan_b]:
            df = (
                pd.DataFrame(
                    dict(
                        pathlength=nx.single_source_dijkstra_path_length(
                            graph, b, cutoff=200
                        )
                    )
                )
                .join(n.buses.country)
                .dropna()
            )
            assert (
                not df.empty
            ), "No buses with defined country within 200km of bus `{}`".format(b)
            n.buses.at[b, "country"] = df.loc[df.pathlength.idxmin(), "country"]

        logger.warning(
            "{} buses are not in any country or offshore shape,"
            " {} have been assigned from the tag of the entsoe map,"
            " the rest from the next bus in terms of pathlength.".format(
                c_nan_b.sum(), c_nan_b.sum() - c_tag_nan_b.sum()
            )
        )

    return buses


def _replace_b2b_converter_at_country_border_by_link(n):
    # Affects only the B2B converter in Lithuania at the Polish border at the moment
    buscntry = n.buses.country
    linkcntry = n.links.bus0.map(buscntry)
    converters_i = n.links.index[
        (n.links.carrier == "B2B") & (linkcntry == n.links.bus1.map(buscntry))
    ]

    def findforeignbus(G, i):
        cntry = linkcntry.at[i]
        for busattr in ("bus0", "bus1"):
            b0 = n.links.at[i, busattr]
            for b1 in G[b0]:
                if buscntry[b1] != cntry:
                    return busattr, b0, b1
        return None, None, None

    for i in converters_i:
        G = n.graph()
        busattr, b0, b1 = findforeignbus(G, i)
        if busattr is not None:
            comp, line = next(iter(G[b0][b1]))
            if comp != "Line":
                logger.warning(
                    "Unable to replace B2B `{}` expected a Line, but found a {}".format(
                        i, comp
                    )
                )
                continue

            n.links.at[i, busattr] = b1
            n.links.at[i, "p_nom"] = min(
                n.links.at[i, "p_nom"], n.lines.at[line, "s_nom"]
            )
            n.links.at[i, "carrier"] = "DC"
            n.links.at[i, "underwater_fraction"] = 0.0
            n.links.at[i, "length"] = n.lines.at[line, "length"]

            n.remove("Line", line)
            n.remove("Bus", b0)

            logger.info(
                "Replacing B2B converter `{}` together with bus `{}` and line `{}` by an HVDC tie-line {}-{}".format(
                    i, b0, line, linkcntry.at[i], buscntry.at[b1]
                )
            )


def _set_links_underwater_fraction(n, offshore_shapes):
    if n.links.empty:
        return

    if not hasattr(n.links, "geometry"):
        n.links["underwater_fraction"] = 0.0
    else:
        offshore_shape = gpd.read_file(offshore_shapes).union_all()
        links = gpd.GeoSeries(n.links.geometry.dropna().map(shapely.wkt.loads))
        n.links["underwater_fraction"] = (
            links.intersection(offshore_shape).length / links.length
        )


def _adjust_capacities_of_under_construction_branches(n, config):
    lines_mode = config["lines"].get("under_construction", "undef")
    if lines_mode == "zero":
        n.lines.loc[n.lines.under_construction, "num_parallel"] = 0.0
        n.lines.loc[n.lines.under_construction, "s_nom"] = 0.0
    elif lines_mode == "remove":
        n.remove("Line", n.lines.index[n.lines.under_construction])
    elif lines_mode != "keep":
        logger.warning(
            "Unrecognized configuration for `lines: under_construction` = `{}`. Keeping under construction lines."
        )

    links_mode = config["links"].get("under_construction", "undef")
    if links_mode == "zero":
        n.links.loc[n.links.under_construction, "p_nom"] = 0.0
    elif links_mode == "remove":
        n.remove("Link", n.links.index[n.links.under_construction])
    elif links_mode != "keep":
        logger.warning(
            "Unrecognized configuration for `links: under_construction` = `{}`. Keeping under construction links."
        )

    if lines_mode == "remove" or links_mode == "remove":
        # We might need to remove further unconnected components
        n = _remove_unconnected_components(n)

    return n


def _set_shapes(n, country_shapes, offshore_shapes):
    # Write the geodataframes country_shapes and offshore_shapes to the network.shapes component
    country_shapes = gpd.read_file(country_shapes).rename(columns={"name": "idx"})
    country_shapes["type"] = "country"
    offshore_shapes = gpd.read_file(offshore_shapes).rename(columns={"name": "idx"})
    offshore_shapes["type"] = "offshore"
    all_shapes = pd.concat([country_shapes, offshore_shapes], ignore_index=True)
    n.add(
        "Shape",
        all_shapes.index,
        geometry=all_shapes.geometry,
        idx=all_shapes.idx,
        type=all_shapes["type"],
    )


def base_network(
    buses,
    converters,
    transformers,
    lines,
    links,
    links_p_nom,
    europe_shape,
    country_shapes,
    offshore_shapes,
    countries,
    parameter_corrections,
    config,
):
    base_network = config["electricity"].get("base_network")
    osm_prebuilt_version = config["electricity"].get("osm-prebuilt-version")
    assert base_network in {
        "entsoegridkit",
        "osm-raw",
        "osm-prebuilt",
    }, f"base_network must be either 'entsoegridkit', 'osm-raw' or 'osm-prebuilt', but got '{base_network}'"
    if base_network == "entsoegridkit":
        warnings.warn(
            "The 'entsoegridkit' base network is deprecated and will be removed in future versions. Please use 'osm-raw' or 'osm-prebuilt' instead.",
            DeprecationWarning,
        )

    logger_str = (
        f"Creating base network using {base_network}"
        + (f" v{osm_prebuilt_version}" if base_network == "osm-prebuilt" else "")
        + "."
    )
    logger.info(logger_str)

    buses = _load_buses(buses, europe_shape, countries, config)
    transformers = _load_transformers(buses, transformers)
    lines = _load_lines(buses, lines)

    if base_network == "entsoegridkit":
        links = _load_links_from_eg(buses, links)
        converters = _load_converters_from_eg(buses, converters)

        # Optionally reconnect Crimea
        if (config["lines"].get("reconnect_crimea", True)) & (
            "UA" in config["countries"]
        ):
            lines = _reconnect_crimea(lines)

        # Set electrical parameters of lines and links
        lines = _set_electrical_parameters_lines_eg(lines, config)
        links = _set_electrical_parameters_links_eg(links, config, links_p_nom)
    elif base_network in {"osm-prebuilt", "osm-raw"}:
        links = _load_links_from_osm(buses, links)
        converters = _load_converters_from_osm(buses, converters)

        # Set electrical parameters of lines and links
        lines = _set_electrical_parameters_lines_osm(lines, config)
        links = _set_electrical_parameters_links_osm(links, config)
    else:
        raise ValueError(
            "base_network must be either 'entsoegridkit', 'osm-raw', or 'osm-prebuilt'"
        )

    # Set electrical parameters of transformers and converters
    transformers = _set_electrical_parameters_transformers(transformers, config)
    converters = _set_electrical_parameters_converters(converters, config)

    n = pypsa.Network()
    n.name = (
        f"PyPSA-Eur ({base_network}"
        + (f" v{osm_prebuilt_version}" if base_network == "osm-prebuilt" else "")
        + ")"
    )

    time = get_snapshots(snakemake.params.snapshots, snakemake.params.drop_leap_day)
    n.set_snapshots(time)

    n.add("Bus", buses.index, **buses)
    n.add("Line", lines.index, **lines)
    n.add("Transformer", transformers.index, **transformers)
    n.add("Link", links.index, **links)
    n.add("Link", converters.index, **converters)

    _set_lines_s_nom_from_linetypes(n)
    if config["electricity"].get("base_network") == "entsoegridkit":
        _apply_parameter_corrections(n, parameter_corrections)

    n = _remove_unconnected_components(n)

    _set_countries_and_substations(n, config, country_shapes, offshore_shapes)

    _set_links_underwater_fraction(n, offshore_shapes)

    _replace_b2b_converter_at_country_border_by_link(n)

    n = _adjust_capacities_of_under_construction_branches(n, config)

    _set_shapes(n, country_shapes, offshore_shapes)

    # Add carriers if they are present in buses.carriers
    carriers_in_buses = set(n.buses.carrier.dropna().unique())
    carriers = carriers_in_buses.intersection({"AC", "DC"})

    if carriers:
        n.add("Carrier", carriers)

    return n


def _get_linetypes_config(line_types, voltages):
    """
    Return the dictionary of linetypes for selected voltages. The dictionary is
    a subset of the dictionary line_types, whose keys match the selected
    voltages.

    Parameters
    ----------
    line_types : dict
        Dictionary of linetypes: keys are nominal voltages and values are linetypes.
    voltages : list
        List of selected voltages.

    Returns
    -------
        Dictionary of linetypes for selected voltages.
    """
    # get voltages value that are not available in the line types
    vnoms_diff = set(voltages).symmetric_difference(set(line_types.keys()))
    if vnoms_diff:
        logger.warning(
            f"Voltages {vnoms_diff} not in the {line_types} or {voltages} list."
        )
    return {k: v for k, v in line_types.items() if k in voltages}


def _get_linetype_by_voltage(v_nom, d_linetypes):
    """
    Return the linetype of a specific line based on its voltage v_nom.

    Parameters
    ----------
    v_nom : float
        The voltage of the line.
    d_linetypes : dict
        Dictionary of linetypes: keys are nominal voltages and values are linetypes.

    Returns
    -------
        The linetype of the line whose nominal voltage is closest to the line voltage.
    """
    v_nom_min, line_type_min = min(
        d_linetypes.items(),
        key=lambda x: abs(x[0] - v_nom),
    )
    return line_type_min


def voronoi(points, outline, crs=4326):
    """
    Create Voronoi polygons from a set of points within an outline.
    """
    pts = gpd.GeoSeries(
        gpd.points_from_xy(points.x, points.y),
        index=points.index,
        crs=crs,
    )
    voronoi = pts.voronoi_polygons(extend_to=outline).clip(outline)

    # can be removed with shapely 2.1 where order is preserved
    # https://github.com/shapely/shapely/issues/2020
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        pts = gpd.GeoDataFrame(geometry=pts)
        voronoi = gpd.GeoDataFrame(geometry=voronoi)
        joined = gpd.sjoin_nearest(pts, voronoi, how="right")

    return joined.dissolve(by="Bus").reindex(points.index).squeeze()


def build_bus_shapes(n, country_shapes, offshore_shapes, countries):
    country_shapes = gpd.read_file(country_shapes).set_index("name")["geometry"]
    offshore_shapes = gpd.read_file(offshore_shapes)
    offshore_shapes = offshore_shapes.reindex(columns=REGION_COLS).set_index("name")[
        "geometry"
    ]

    onshore_regions = []
    offshore_regions = []

    for country in countries:
        c_b = n.buses.country == country

        onshore_shape = country_shapes[country]
        onshore_locs = (
            n.buses.loc[c_b & n.buses.onshore_bus]
            .sort_values(
                by="substation_lv", ascending=False
            )  # preference for substations
            .drop_duplicates(subset=["x", "y"], keep="first")[["x", "y"]]
        )
        onshore_regions.append(
            gpd.GeoDataFrame(
                {
                    "name": onshore_locs.index,
                    "x": onshore_locs["x"],
                    "y": onshore_locs["y"],
                    "geometry": voronoi(onshore_locs, onshore_shape),
                    "country": country,
                },
                crs=n.crs,
            )
        )

        if country not in offshore_shapes.index:
            continue
        offshore_shape = offshore_shapes[country]
        offshore_locs = n.buses.loc[c_b & n.buses.substation_off, ["x", "y"]]
        offshore_regions_c = gpd.GeoDataFrame(
            {
                "name": offshore_locs.index,
                "x": offshore_locs["x"],
                "y": offshore_locs["y"],
                "geometry": voronoi(offshore_locs, offshore_shape),
                "country": country,
            },
            crs=n.crs,
        )
        sel = offshore_regions_c.to_crs(3035).area > 10  # m2
        offshore_regions_c = offshore_regions_c.loc[sel]
        offshore_regions.append(offshore_regions_c)

    shapes = pd.concat(onshore_regions, ignore_index=True).set_crs(n.crs)

    return onshore_regions, offshore_regions, shapes, offshore_shapes


def append_bus_shapes(n, shapes, type):
    """
    Append shapes to the network. If shapes with the same component and type
    already exist, they will be removed.

    Parameters:
        n (pypsa.Network): The network to which the shapes will be appended.
        shapes (geopandas.GeoDataFrame): The shapes to be appended.
        **kwargs: Additional keyword arguments used in `n.add`.

    Returns:
        None
    """
    remove = n.shapes.query("component == 'Bus' and type == @type").index
    n.remove("Shape", remove)

    offset = n.shapes.index.astype(int).max() + 1 if not n.shapes.empty else 0
    shapes = shapes.rename(lambda x: int(x) + offset)
    n.add(
        "Shape",
        shapes.index,
        geometry=shapes.geometry,
        idx=shapes.name,
        component="Bus",
        type=type,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("base_network")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    countries = snakemake.params.countries

    buses = snakemake.input.buses
    converters = snakemake.input.converters
    transformers = snakemake.input.transformers
    lines = snakemake.input.lines
    links = snakemake.input.links
    europe_shape = snakemake.input.europe_shape
    country_shapes = snakemake.input.country_shapes
    offshore_shapes = snakemake.input.offshore_shapes
    config = snakemake.config

    if "links_p_nom" in snakemake.input.keys():
        links_p_nom = snakemake.input.links_p_nom
    else:
        links_p_nom = None

    if "parameter_corrections" in snakemake.input.keys():
        parameter_corrections = snakemake.input.parameter_corrections
    else:
        parameter_corrections = None

    n = base_network(
        buses,
        converters,
        transformers,
        lines,
        links,
        links_p_nom,
        europe_shape,
        country_shapes,
        offshore_shapes,
        countries,
        parameter_corrections,
        config,
    )

    onshore_regions, offshore_regions, shapes, offshore_shapes = build_bus_shapes(
        n,
        country_shapes,
        offshore_shapes,
        countries,
    )

    shapes.to_file(snakemake.output.regions_onshore)
    # append_bus_shapes(n, shapes, "onshore")

    if offshore_regions:
        shapes = pd.concat(offshore_regions, ignore_index=True)
        shapes.to_file(snakemake.output.regions_offshore)
        # append_bus_shapes(n, shapes, "offshore")
    else:
        offshore_shapes.to_frame().to_file(snakemake.output.regions_offshore)

    n.meta = snakemake.config
    n.export_to_netcdf(snakemake.output.base_network)
