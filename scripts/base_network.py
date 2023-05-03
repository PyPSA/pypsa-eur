# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8
"""
Creates the network topology from a `ENTSO-E map extract.

<https://github.com/PyPSA/GridKit/tree/master/entsoe>`_ (March 2022) as a PyPSA
network.

Relevant Settings
-----------------

.. code:: yaml

    snapshots:

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
        include_tyndp:

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
- ``data/links_tyndp.csv``: List of projects in the `TYNDP 2018 <https://tyndp.entsoe.eu/tyndp2018/>`_ that are at least *in permitting* with fields for start- and endpoint (names and coordinates), length, capacity, construction status, and project reference ID.
- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``resources/offshore_shapes.geojson``: confer :ref:`shapes`
- ``resources/europe_shape.geojson``: confer :ref:`shapes`

Outputs
-------

- ``networks/base.nc``

    .. image:: img/base.png
        :scale: 33 %

Description
-----------
"""

import logging
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
from _helpers import configure_logging
from scipy import spatial
from scipy.sparse import csgraph
from shapely.geometry import LineString, Point

logger = logging.getLogger(__name__)


def _get_oid(df):
    if "tags" in df.columns:
        return df.tags.str.extract('"oid"=>"(\d+)"', expand=False)
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
    tree = spatial.KDTree(treecoords)
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


def _load_buses_from_eg(eg_buses, europe_shape, config_elec):
    buses = (
        pd.read_csv(
            eg_buses,
            quotechar="'",
            true_values=["t"],
            false_values=["f"],
            dtype=dict(bus_id="str"),
        )
        .set_index("bus_id")
        .drop(["station_id"], axis=1)
        .rename(columns=dict(voltage="v_nom"))
    )

    buses["carrier"] = buses.pop("dc").map({True: "DC", False: "AC"})
    buses["under_construction"] = buses["under_construction"].fillna(False).astype(bool)

    # remove all buses outside of all countries including exclusive economic zones (offshore)
    europe_shape = gpd.read_file(europe_shape).loc[0, "geometry"]
    europe_shape_prepped = shapely.prepared.prep(europe_shape)
    buses_in_europe_b = buses[["x", "y"]].apply(
        lambda p: europe_shape_prepped.contains(Point(p)), axis=1
    )

    buses_with_v_nom_to_keep_b = (
        buses.v_nom.isin(config_elec["voltages"]) | buses.v_nom.isnull()
    )
    logger.info(
        "Removing buses with voltages {}".format(
            pd.Index(buses.v_nom.unique()).dropna().difference(config_elec["voltages"])
        )
    )

    return pd.DataFrame(buses.loc[buses_in_europe_b & buses_with_v_nom_to_keep_b])


def _load_transformers_from_eg(buses, eg_transformers):
    transformers = pd.read_csv(
        eg_transformers,
        quotechar="'",
        true_values=["t"],
        false_values=["f"],
        dtype=dict(transformer_id="str", bus0="str", bus1="str"),
    ).set_index("transformer_id")

    transformers = _remove_dangling_branches(transformers, buses)

    return transformers


def _load_converters_from_eg(buses, eg_converters):
    converters = pd.read_csv(
        eg_converters,
        quotechar="'",
        true_values=["t"],
        false_values=["f"],
        dtype=dict(converter_id="str", bus0="str", bus1="str"),
    ).set_index("converter_id")

    converters = _remove_dangling_branches(converters, buses)

    converters["carrier"] = "B2B"

    return converters


def _load_links_from_eg(buses, eg_links):
    links = pd.read_csv(
        eg_links,
        quotechar="'",
        true_values=["t"],
        false_values=["f"],
        dtype=dict(link_id="str", bus0="str", bus1="str", under_construction="bool"),
    ).set_index("link_id")

    links["length"] /= 1e3

    # Skagerrak Link is connected to 132kV bus which is removed in _load_buses_from_eg.
    # Connect to neighboring 380kV bus
    links.loc[links.bus1 == "6396", "bus1"] = "6398"

    links = _remove_dangling_branches(links, buses)

    # Add DC line parameters
    links["carrier"] = "DC"

    return links


def _add_links_from_tyndp(buses, links, links_tyndp, europe_shape):
    links_tyndp = pd.read_csv(links_tyndp)

    # remove all links from list which lie outside all of the desired countries
    europe_shape = gpd.read_file(europe_shape).loc[0, "geometry"]
    europe_shape_prepped = shapely.prepared.prep(europe_shape)
    x1y1_in_europe_b = links_tyndp[["x1", "y1"]].apply(
        lambda p: europe_shape_prepped.contains(Point(p)), axis=1
    )
    x2y2_in_europe_b = links_tyndp[["x2", "y2"]].apply(
        lambda p: europe_shape_prepped.contains(Point(p)), axis=1
    )
    is_within_covered_countries_b = x1y1_in_europe_b & x2y2_in_europe_b

    if not is_within_covered_countries_b.all():
        logger.info(
            "TYNDP links outside of the covered area (skipping): "
            + ", ".join(links_tyndp.loc[~is_within_covered_countries_b, "Name"])
        )

        links_tyndp = links_tyndp.loc[is_within_covered_countries_b]
        if links_tyndp.empty:
            return buses, links

    has_replaces_b = links_tyndp.replaces.notnull()
    oids = dict(Bus=_get_oid(buses), Link=_get_oid(links))
    keep_b = dict(
        Bus=pd.Series(True, index=buses.index), Link=pd.Series(True, index=links.index)
    )
    for reps in links_tyndp.loc[has_replaces_b, "replaces"]:
        for comps in reps.split(":"):
            oids_to_remove = comps.split(".")
            c = oids_to_remove.pop(0)
            keep_b[c] &= ~oids[c].isin(oids_to_remove)
    buses = buses.loc[keep_b["Bus"]]
    links = links.loc[keep_b["Link"]]

    links_tyndp["j"] = _find_closest_links(
        links, links_tyndp, distance_upper_bound=0.20
    )
    # Corresponds approximately to 20km tolerances

    if links_tyndp["j"].notnull().any():
        logger.info(
            "TYNDP links already in the dataset (skipping): "
            + ", ".join(links_tyndp.loc[links_tyndp["j"].notnull(), "Name"])
        )
        links_tyndp = links_tyndp.loc[links_tyndp["j"].isnull()]
        if links_tyndp.empty:
            return buses, links

    tree = spatial.KDTree(buses[["x", "y"]])
    _, ind0 = tree.query(links_tyndp[["x1", "y1"]])
    ind0_b = ind0 < len(buses)
    links_tyndp.loc[ind0_b, "bus0"] = buses.index[ind0[ind0_b]]

    _, ind1 = tree.query(links_tyndp[["x2", "y2"]])
    ind1_b = ind1 < len(buses)
    links_tyndp.loc[ind1_b, "bus1"] = buses.index[ind1[ind1_b]]

    links_tyndp_located_b = (
        links_tyndp["bus0"].notnull() & links_tyndp["bus1"].notnull()
    )
    if not links_tyndp_located_b.all():
        logger.warning(
            "Did not find connected buses for TYNDP links (skipping): "
            + ", ".join(links_tyndp.loc[~links_tyndp_located_b, "Name"])
        )
        links_tyndp = links_tyndp.loc[links_tyndp_located_b]

    logger.info("Adding the following TYNDP links: " + ", ".join(links_tyndp["Name"]))

    links_tyndp = links_tyndp[["bus0", "bus1"]].assign(
        carrier="DC",
        p_nom=links_tyndp["Power (MW)"],
        length=links_tyndp["Length (given) (km)"].fillna(
            links_tyndp["Length (distance*1.2) (km)"]
        ),
        under_construction=True,
        underground=False,
        geometry=(
            links_tyndp[["x1", "y1", "x2", "y2"]].apply(
                lambda s: str(LineString([[s.x1, s.y1], [s.x2, s.y2]])), axis=1
            )
        ),
        tags=(
            '"name"=>"'
            + links_tyndp["Name"]
            + '", '
            + '"ref"=>"'
            + links_tyndp["Ref"]
            + '", '
            + '"status"=>"'
            + links_tyndp["status"]
            + '"'
        ),
    )

    links_tyndp.index = "T" + links_tyndp.index.astype(str)

    links = pd.concat([links, links_tyndp], sort=True)

    return buses, links


def _load_lines_from_eg(buses, eg_lines):
    lines = (
        pd.read_csv(
            eg_lines,
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


def _set_electrical_parameters_lines(lines, config):
    v_noms = config["electricity"]["voltages"]
    linetypes = config["lines"]["types"]

    for v_nom in v_noms:
        lines.loc[lines["v_nom"] == v_nom, "type"] = linetypes[v_nom]

    lines["s_max_pu"] = config["lines"]["s_max_pu"]

    return lines


def _set_lines_s_nom_from_linetypes(n):
    n.lines["s_nom"] = (
        np.sqrt(3)
        * n.lines["type"].map(n.line_types.i_nom)
        * n.lines["v_nom"]
        * n.lines.num_parallel
    )


def _set_electrical_parameters_links(links, config, links_p_nom):
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


def _set_electrical_parameters_converters(converters, config):
    p_max_pu = config["links"].get("p_max_pu", 1.0)
    converters["p_max_pu"] = p_max_pu
    converters["p_min_pu"] = -p_max_pu

    converters["p_nom"] = 2000

    # Converters are combined with links
    converters["under_construction"] = False
    converters["underground"] = False

    return converters


def _set_electrical_parameters_transformers(transformers, config):
    config = config["transformers"]

    ## Add transformer parameters
    transformers["x"] = config.get("x", 0.1)
    transformers["s_nom"] = config.get("s_nom", 2000)
    transformers["type"] = config.get("type", "")

    return transformers


def _remove_dangling_branches(branches, buses):
    return pd.DataFrame(
        branches.loc[branches.bus0.isin(buses.index) & branches.bus1.isin(buses.index)]
    )


def _remove_unconnected_components(network):
    _, labels = csgraph.connected_components(network.adjacency_matrix(), directed=False)
    component = pd.Series(labels, index=network.buses.index)

    component_sizes = component.value_counts()
    components_to_remove = component_sizes.iloc[1:]

    logger.info(
        "Removing {} unconnected network components with less than {} buses. In total {} buses.".format(
            len(components_to_remove),
            components_to_remove.max(),
            components_to_remove.sum(),
        )
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

    gb = buses.loc[substation_b].groupby(
        ["x", "y"], as_index=False, group_keys=False, sort=False
    )
    bus_map_low = gb.apply(prefer_voltage, "min")
    lv_b = (bus_map_low == bus_map_low.index).reindex(buses.index, fill_value=False)
    bus_map_high = gb.apply(prefer_voltage, "max")
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

    buses["substation_lv"] = (
        lv_b & onshore_b & (~buses["under_construction"]) & has_connections_b
    )
    buses["substation_off"] = (offshore_b | (hv_b & onshore_b)) & (
        ~buses["under_construction"]
    )

    c_nan_b = buses.country.isnull()
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
        offshore_shape = gpd.read_file(offshore_shapes).unary_union
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
        n.mremove("Line", n.lines.index[n.lines.under_construction])
    elif lines_mode != "keep":
        logger.warning(
            "Unrecognized configuration for `lines: under_construction` = `{}`. Keeping under construction lines."
        )

    links_mode = config["links"].get("under_construction", "undef")
    if links_mode == "zero":
        n.links.loc[n.links.under_construction, "p_nom"] = 0.0
    elif links_mode == "remove":
        n.mremove("Link", n.links.index[n.links.under_construction])
    elif links_mode != "keep":
        logger.warning(
            "Unrecognized configuration for `links: under_construction` = `{}`. Keeping under construction links."
        )

    if lines_mode == "remove" or links_mode == "remove":
        # We might need to remove further unconnected components
        n = _remove_unconnected_components(n)

    return n


def base_network(
    eg_buses,
    eg_converters,
    eg_transformers,
    eg_lines,
    eg_links,
    links_p_nom,
    links_tyndp,
    europe_shape,
    country_shapes,
    offshore_shapes,
    parameter_corrections,
    config,
):
    buses = _load_buses_from_eg(eg_buses, europe_shape, config["electricity"])

    links = _load_links_from_eg(buses, eg_links)
    if config["links"].get("include_tyndp"):
        buses, links = _add_links_from_tyndp(buses, links, links_tyndp, europe_shape)

    converters = _load_converters_from_eg(buses, eg_converters)

    lines = _load_lines_from_eg(buses, eg_lines)
    transformers = _load_transformers_from_eg(buses, eg_transformers)

    lines = _set_electrical_parameters_lines(lines, config)
    transformers = _set_electrical_parameters_transformers(transformers, config)
    links = _set_electrical_parameters_links(links, config, links_p_nom)
    converters = _set_electrical_parameters_converters(converters, config)

    n = pypsa.Network()
    n.name = "PyPSA-Eur"

    n.set_snapshots(pd.date_range(freq="h", **config["snapshots"]))

    n.import_components_from_dataframe(buses, "Bus")
    n.import_components_from_dataframe(lines, "Line")
    n.import_components_from_dataframe(transformers, "Transformer")
    n.import_components_from_dataframe(links, "Link")
    n.import_components_from_dataframe(converters, "Link")

    _set_lines_s_nom_from_linetypes(n)

    _apply_parameter_corrections(n, parameter_corrections)

    n = _remove_unconnected_components(n)

    _set_countries_and_substations(n, config, country_shapes, offshore_shapes)

    _set_links_underwater_fraction(n, offshore_shapes)

    _replace_b2b_converter_at_country_border_by_link(n)

    n = _adjust_capacities_of_under_construction_branches(n, config)

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("base_network")
    configure_logging(snakemake)

    n = base_network(
        snakemake.input.eg_buses,
        snakemake.input.eg_converters,
        snakemake.input.eg_transformers,
        snakemake.input.eg_lines,
        snakemake.input.eg_links,
        snakemake.input.links_p_nom,
        snakemake.input.links_tyndp,
        snakemake.input.europe_shape,
        snakemake.input.country_shapes,
        snakemake.input.offshore_shapes,
        snakemake.input.parameter_corrections,
        snakemake.config,
    )

    n.meta = snakemake.config
    n.export_to_netcdf(snakemake.output[0])
