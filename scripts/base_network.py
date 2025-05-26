# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


"""
Creates the network topology from a `ENTSO-E map extract.
<https://github.com/PyPSA/GridKit/tree/master/entsoe>`_ (March 2022)
or `OpenStreetMap data <https://www.openstreetmap.org/>`_ (Aug 2024)
as a PyPSA
network.

Description
-----------
Creates the network topology from an ENTSO-E map extract, and create Voronoi shapes for each bus representing both onshore and offshore regions.
"""

import logging
import multiprocessing as mp
import warnings
from functools import partial
from itertools import chain, product

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pypsa
import shapely
import shapely.prepared
import shapely.wkt
import yaml
from packaging.version import Version, parse
from scipy.sparse import csgraph
from scipy.spatial import KDTree
from shapely.geometry import Point
from tqdm import tqdm

from scripts._helpers import (
    REGION_COLS,
    configure_logging,
    get_snapshots,
    set_scenario_config,
)

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


def _load_converters_from_raw(buses, converters):
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
    # Connect to neighbouring 380kV bus
    links.loc[links.bus1 == "6396", "bus1"] = "6398"

    links = _remove_dangling_branches(links, buses)

    # Add DC line parameters
    links["carrier"] = "DC"

    return links


def _load_links_from_raw(buses, links):
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


def _set_electrical_parameters_lines_raw(lines, config):
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
    p_min_pu = config["links"].get("p_min_pu", -p_max_pu)
    links["p_max_pu"] = p_max_pu
    links["p_min_pu"] = p_min_pu

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


def _set_electrical_parameters_links_raw(links, config):
    if links.empty:
        return links

    p_max_pu = config["links"].get("p_max_pu", 1.0)
    p_min_pu = config["links"].get("p_min_pu", -p_max_pu)
    links["p_max_pu"] = p_max_pu
    links["p_min_pu"] = p_min_pu
    links["carrier"] = "DC"
    links["dc"] = True

    return links


def _set_electrical_parameters_converters(converters, config):
    p_max_pu = config["links"].get("p_max_pu", 1.0)
    p_min_pu = config["links"].get("p_min_pu", -p_max_pu)
    converters["p_max_pu"] = p_max_pu
    converters["p_min_pu"] = p_min_pu

    # if column "p_nom" does not exist, set to 2000
    if "p_nom" not in converters:
        converters["p_nom"] = 2000

    # Converters are combined with links
    converters["under_construction"] = False
    converters["underground"] = False
    converters["dc"] = False  # ToDo Find a better assumption

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
            assert not df.empty, (
                f"No buses with defined country within 200km of bus `{b}`"
            )
            n.buses.at[b, "country"] = df.loc[df.pathlength.idxmin(), "country"]

        logger.warning(
            f"{c_nan_b.sum()} buses are not in any country or offshore shape,"
            f" {c_nan_b.sum() - c_tag_nan_b.sum()} have been assigned from the tag of the entsoe map,"
            " the rest from the next bus in terms of pathlength."
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
                    f"Unable to replace B2B `{i}` expected a Line, but found a {comp}"
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
                f"Replacing B2B converter `{i}` together with bus `{b0}` and line `{line}` by an HVDC tie-line {linkcntry.at[i]}-{buscntry.at[b1]}"
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
    assert base_network in {"entsoegridkit", "osm-raw", "osm-prebuilt", "tyndp-raw"}, (
        f"base_network must be either 'entsoegridkit', 'osm-raw', 'osm-prebuilt' or 'tyndp-raw', but got '{base_network}'"
    )
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
    elif base_network in {"osm-prebuilt", "osm-raw", "tyndp-raw"}:
        links = _load_links_from_raw(buses, links)
        converters = _load_converters_from_raw(buses, converters)

        # Set electrical parameters of lines and links
        lines = _set_electrical_parameters_lines_raw(lines, config)
        links = _set_electrical_parameters_links_raw(links, config)
    else:
        raise ValueError(
            "base_network must be either 'entsoegridkit', 'osm-raw', 'osm-prebuilt', or 'tyndp-raw'"
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


def process_onshore_regions(
    adm: str,
    buses: pd.DataFrame,
    admin_shapes: gpd.GeoDataFrame,
    crs: str,
) -> gpd.GeoDataFrame:
    country = admin_shapes.loc[adm, "country"]
    c_b = buses.admin == adm

    onshore_shape = admin_shapes.loc[adm, "geometry"]
    onshore_locs = (
        buses.loc[c_b & buses.onshore_bus]
        .sort_values(by="substation_lv", ascending=False)  # preference for substations
        .drop_duplicates(subset=["x", "y", "country"], keep="first")[
            ["x", "y", "country"]
        ]
    )
    onshore_regions_adm = gpd.GeoDataFrame(
        {
            "name": onshore_locs.index,
            "x": onshore_locs["x"],
            "y": onshore_locs["y"],
            "geometry": voronoi(onshore_locs, onshore_shape),
            "country": country,
        },
        crs=crs,
    )

    return onshore_regions_adm


def process_offshore_regions(
    buses: pd.DataFrame,
    offshore_shapes: gpd.GeoDataFrame,
    countries: list[str],
    crs: str,
) -> list[gpd.GeoDataFrame]:
    offshore_regions = []

    tqdm_kwargs = dict(
        ascii=False,
        unit=" regions",
        total=len(countries),
        desc="Building offshore regions",
    )
    for country in tqdm(countries, **tqdm_kwargs):
        if country not in offshore_shapes.index:
            continue

        c_b = buses.country == country
        offshore_shape = offshore_shapes[country]
        offshore_locs = buses.loc[c_b & buses.substation_off, ["x", "y"]]
        offshore_regions_c = gpd.GeoDataFrame(
            {
                "name": offshore_locs.index,
                "x": offshore_locs["x"],
                "y": offshore_locs["y"],
                "geometry": voronoi(offshore_locs, offshore_shape),
                "country": country,
            },
            crs=crs,
        )
        sel = offshore_regions_c.to_crs(3035).area > 10  # m2
        offshore_regions_c = offshore_regions_c.loc[sel]
        offshore_regions.append(offshore_regions_c)

    return offshore_regions


def build_bus_shapes(
    n: pypsa.Network,
    admin_shapes: gpd.GeoDataFrame,
    offshore_shapes: str,
    countries: list[str],
) -> tuple[
    list[gpd.GeoDataFrame], list[gpd.GeoDataFrame], gpd.GeoDataFrame, gpd.GeoDataFrame
]:
    """
    Build onshore and offshore regions for buses in the network.

    Parameters
    ----------
        n (pypsa.Network) : The network for which the bus shapes will be built.
        admin_shapes (gpd.GeoDataFrame) : GeoDataFrame with administrative region shapes indexed by name.
        offshore_shapes (str) : Path to the file containing offshore shapes.
        countries (list[str]) : List of country codes to process.

    Returns
    -------
        tuple[list[gpd.GeoDataFrame], list[gpd.GeoDataFrame], gpd.GeoDataFrame, gpd.GeoDataFrame]

        A tuple containing:
            - List of GeoDataFrames for each onshore region
            - List of GeoDataFrames for each offshore region
            - Combined GeoDataFrame of all onshore shapes
            - Combined GeoDataFrame of all offshore shapes
    """
    offshore_shapes = gpd.read_file(offshore_shapes)
    offshore_shapes = offshore_shapes.reindex(columns=REGION_COLS).set_index("name")[
        "geometry"
    ]

    buses = n.buses[
        ["x", "y", "country", "onshore_bus", "substation_lv", "substation_off"]
    ].copy()

    buses["geometry"] = gpd.points_from_xy(buses["x"], buses["y"])
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs="EPSG:4326")
    buses["admin"] = ""

    # Map buses per country
    for country in countries:
        buses_subset = buses.loc[buses["country"] == country]

        buses.loc[buses_subset.index, "admin"] = gpd.sjoin_nearest(
            buses_subset.to_crs(epsg=3857),
            admin_shapes.loc[admin_shapes["country"] == country].to_crs(epsg=3857),
            how="left",
        )["admin_right"]

    # Create Voronoi polygons for each administrative region.
    # If administrative clustering is deactivated, voronoi cells are created on a country level.
    admin_regions = sorted(
        set(buses.admin.unique()).intersection(admin_shapes.index.unique())
    )

    # Onshore regions
    nprocesses = snakemake.threads
    tqdm_kwargs = dict(
        ascii=False,
        unit=" regions",
        total=len(admin_regions),
        desc="Building onshore regions",
    )
    func = partial(
        process_onshore_regions,
        buses=buses,
        admin_shapes=admin_shapes,
        crs=n.crs.name,
    )

    with mp.Pool(processes=nprocesses) as pool:
        onshore_regions = list(tqdm(pool.imap(func, admin_regions), **tqdm_kwargs))
    onshore_shapes = pd.concat(onshore_regions, ignore_index=True).set_crs(n.crs)
    logger.info(f"In total {len(onshore_shapes)} onshore regions.")

    # Offshore regions
    offshore_regions = process_offshore_regions(
        buses,
        offshore_shapes,
        countries,
        n.crs.name,
    )
    if offshore_regions:
        offshore_shapes = pd.concat(offshore_regions, ignore_index=True).set_crs(n.crs)
    else:
        offshore_shapes = gpd.GeoDataFrame(
            columns=["name", "geometry"], crs=n.crs
        ).set_index("name")
    logger.info(f"In total {len(offshore_shapes)} offshore regions.")

    return onshore_regions, offshore_regions, onshore_shapes, offshore_shapes


def append_bus_shapes(n, shapes, type):
    """
    Append shapes to the network. If shapes with the same component and type
    already exist, they will be removed.

    Parameters
    ----------
        n (pypsa.Network): The network to which the shapes will be appended.
        shapes (geopandas.GeoDataFrame): The shapes to be appended.
        **kwargs: Additional keyword arguments used in `n.add`.

    Returns
    -------
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


def find_neighbours(
    polygon: shapely.geometry.Polygon, index: str, gdf: gpd.GeoDataFrame
) -> list:
    """
    Find neighbouring polygons of a given polygon in a GeoDataFrame using spatial index.

    Parameters
    ----------
        polygon (shapely.geometry.Polygon): Polygon for which to find neighbours.
        index (str): Index of the polygon.
        gdf (gpd.GeoDataFrame): GeoDataFrame containing all polygons.

    Returns
    -------
        list: List of indices of neighbouring polygons.
    """
    possible_neighbours = gdf.sindex.intersection(polygon.bounds)

    # Get actual neighbours by filtering those that touch
    neighbours = gdf.iloc[list(possible_neighbours)]
    neighbours = neighbours[neighbours.geometry.touches(polygon)]

    # if no neighbours exist, return empty list
    if len(neighbours) == 0:
        return None

    # Return the indices of touching neighbours
    return neighbours.index.tolist()


def keep_good_neighbours(
    adm: str,
    neighbours: list,
    parent_dict: dict,
    country_dict: dict,
) -> list:
    """
    Filtering list of neighbours based on country and parent.

    Parameters
    ----------
        adm (str): Index of the administrative region.
        neighbours (list): List of neighbours.
        parent_dict (dict): Dictionary with parent of each administrative region.
        country_dict (dict): Dictionary with country of each administrative region.

    Returns
    -------
        list: List of filtered
    """
    # Only keep neighbours that are located in the same country

    neighbours = [n for n in neighbours if n[:2] == country_dict[adm]]

    # Filter new neighbours by parent
    new_neighbours = [n for n in neighbours if parent_dict[n] == parent_dict[adm]]

    # If no neighbours are left, keep all neighbours
    if not new_neighbours:
        return neighbours

    return new_neighbours


def sort_values_by_dict(
    neighbours: list,
    dicts: list,
    ascending: bool = True,
) -> list:
    """
    Sorts a list of keys by values from multiple dictionaries in order of priority.

    Parameters
    ----------
        neighbours (list): List of keys to sort.
        dicts (list): List of dictionaries containing values to sort by.
        ascending (bool): Whether to sort in ascending order.

    Returns
    -------
        list: Sorted list of keys.
    """
    return sorted(
        neighbours,
        key=lambda x: tuple(
            d[x] for d in dicts
        ),  # Create sorting tuple from multiple dictionaries
        reverse=not ascending,
    )


def create_merged_admin_region(
    row: pd.Series,
    first_neighbours_dict: dict,
    admin_shapes: gpd.GeoDataFrame,
) -> pd.Series:
    """
    Creates a merged administrative region from a given row and its first neighbours.

    Parameters
    ----------
        row (pd.Series): Series containing information about the region to be merged.
        first_neighbours_dict (dict): Dictionary containing first neighbours for each region.
        admin_shapes (gpd.GeoDataFrame): GeoDataFrame containing all administrative regions.

    Returns
    -------
        pd.Series: Series containing information about the merged region.
    """
    first_neighbours = first_neighbours_dict[row.name]
    neighbours_contain = list(
        set(
            chain.from_iterable(
                [
                    admin_shapes.loc[neighbour, "contains"]
                    for neighbour in first_neighbours
                ]
            )
        )
    )

    merge_regions = sorted([row.name] + first_neighbours)
    contains_cumulative = sorted(
        set(row.contains + first_neighbours + neighbours_contain)
    )
    gdf = admin_shapes.loc[merge_regions]

    geometry = gdf.union_all()
    country = gdf.loc[row.name, "country"]
    parent = gdf.loc[row.name, "parent"]
    substations = gdf["substations"].sum()
    isempty = not any(
        ~gdf["isempty"].values
    )  # if one of the regions is not empty, the merged region is not empty
    area = gdf["area"].sum()
    neighbours = sorted(
        list(
            set(
                chain.from_iterable(
                    [
                        neighbour
                        for neighbour in gdf["neighbours"]
                        if isinstance(neighbour, list)
                    ]
                )
            )
        )
    )
    neighbours = [
        n for n in neighbours if n not in merge_regions
    ]  # Remove contained values from neighbours

    return pd.Series(
        {
            "geometry": geometry,
            "country": country,
            "parent": parent,
            "substations": substations,
            "contains": contains_cumulative,  # List of contains
            "isempty": isempty,
            "area": area,
            "neighbours": neighbours,  # List of neighbours
        }
    )


def update_names(
    names: list[str],
) -> str:
    """
    Updating names of merged administrative regions using the alphabetically first name in a list.

    Parameters
    ----------
        names (list): List of names to update.

    Returns
    -------
        str: Updated name.
    """
    if len(names) == 1:
        return names[0]

    basename = sorted(names)[0]
    name = basename + "+" + str(len(names) - 1)

    return name


def clean_dict(
    diction: dict,
) -> dict:
    """
    Cleans a dictionary by removing duplicates (keys or values that appear in multiple key-value pairs).

    Parameters
    ----------
        diction (dict): Dictionary to clean.

    Returns
    -------
        dict: Cleaned dictionary.
    """

    if not diction:
        return diction

    # Sub dictionary where values are only of length 1
    single_occurrences = {k: v for k, v in diction.items() if len(v) == 1}

    # Enters only if single_occurrences is not empty
    if single_occurrences:
        # Create list of key, value pairs and store in list
        single_values = list(chain.from_iterable(single_occurrences.values()))
        single_keys = list(single_occurrences.keys())

        # Create dict with single_keys() as keys and sets of single_keys and single_values as values using
        single_df = pd.DataFrame(
            zip(
                single_keys,
                [set([key, value]) for key, value in zip(single_keys, single_values)],
            )
        )
        # Only keep first occurrence of duplicates
        single_df = single_df.drop_duplicates(subset=1, keep="first")

        # Difference between single_df[0].values and single_keys
        drop_keys = set(single_keys) - set(single_df[0].values)

        # Drop keys in double_values
        diction = {k: v for k, v in diction.items() if k not in drop_keys}

    # Create list of values, that also exist as key
    double_values = list(chain.from_iterable(diction.values()))
    double_values = [x for x in double_values if x in diction.keys()]

    for key, values in diction.items():
        diction[key] = [value for value in values if value not in double_values]

    # Remove keys that have no values left
    diction = {k: v for k, v in diction.items() if v}

    return diction


def get_nearest_neighbour(
    row: pd.Series,
    admin_shapes: gpd.GeoDataFrame,
) -> str:
    """
    Finds the nearest neighbour containing substations within the same region.

    Parameters
    ----------
        row (pd.Series): Series containing information about the region.
        admin_shapes (gpd.GeoDataFrame): GeoDataFrame containing all administrative regions.

    Returns
    -------
        str: Index of the nearest neighbour.
    """
    country = row["country"]
    gdf = gpd.GeoDataFrame([row.loc[["country", "geometry"]]], crs=admin_shapes.crs)
    nearest_neighbours = admin_shapes.loc[
        (admin_shapes.index != row.name)
        & (admin_shapes["country"] == country)
        & (admin_shapes["isempty"] == False),
        ["country", "geometry"],
    ]

    try:
        nearest_neighbour = (
            gdf.to_crs(epsg=3035)
            .sjoin_nearest(
                nearest_neighbours.to_crs(epsg=3035),
                how="left",
            )["index_right"]
            .values[0]
        )
    except KeyError:
        logger.warning(f"Skipping for {country}")
        nearest_neighbour = ""

    return nearest_neighbour


def merge_regions_recursive(
    admin_shapes: gpd.GeoDataFrame,
    neighbours_missing: bool = True,
) -> gpd.GeoDataFrame:
    """
    Recursive function to merge administrative regions that do not contain substations with their neighbours.
    Prioritises neighbours with the most substations and smallest area.
    Terminates when all regions without substations have been merged.

    Parameters
    ----------
        admin_shapes (gpd.GeoDataFrame): GeoDataFrame containing all administrative regions.
        neighbours_missing (bool): Whether to find neighbours if they are missing.

    Returns
    -------
        gpd.GeoDataFrame: GeoDataFrame containing the merged administrative regions
    """
    while True:
        # Calculate area
        admin_shapes["area"] = admin_shapes.to_crs(epsg=3035).area
        area_dict = admin_shapes["area"].to_dict()
        country_dict = admin_shapes["country"].to_dict()
        parent_dict = admin_shapes["parent"].to_dict()
        substations_dict = (
            admin_shapes["substations"].mul(-1).to_dict()
        )  # multiply by -1 to sort in ascending order

        if neighbours_missing:
            # Find all neighbours
            admin_shapes["neighbours"] = None
            admin_shapes.loc[admin_shapes.isempty, "neighbours"] = admin_shapes.loc[
                admin_shapes.isempty
            ].apply(
                lambda row: find_neighbours(row.geometry, row.name, admin_shapes),
                axis=1,
            )

            # Keep only viable neighbours (preferably same parent)
            admin_shapes.loc[
                admin_shapes.isempty & ~admin_shapes.neighbours.isna(), "neighbours"
            ] = admin_shapes.loc[
                admin_shapes.isempty & ~admin_shapes.neighbours.isna()
            ].apply(
                lambda row: keep_good_neighbours(
                    row.name, row.neighbours, parent_dict, country_dict
                ),
                axis=1,
            )

        # Filter for regions that have neighbours and are empty
        b_isempty_hasneighbours = admin_shapes.isempty & admin_shapes.neighbours.apply(
            lambda x: len(x) > 0 if isinstance(x, list) else False
        )

        if not b_isempty_hasneighbours.any():
            logger.info("All administrative regions without buses have been merged.")
            break

        # Sort neighbours in such a way, that the neighbour with the most substations is first.
        # If there are multiple neighbours with the same number of substations, the one with the smallest area is first.
        # ascending = True (default) - For this to work, substations_dict was multiplied by -1
        # Most negative number of substation (smallest) comes first
        admin_shapes.loc[b_isempty_hasneighbours, "neighbours"] = admin_shapes.loc[
            b_isempty_hasneighbours
        ].apply(
            lambda row: sort_values_by_dict(
                row.neighbours, [substations_dict, area_dict]
            ),
            axis=1,
        )

        # Find all first neighbours and
        # create dict with first neighbours as keys and list of regions of which they are neighbours as values
        first_neighbours = admin_shapes.loc[
            b_isempty_hasneighbours, "neighbours"
        ].apply(lambda x: x[0])
        first_neighbours = (
            first_neighbours.groupby(first_neighbours)
            .apply(lambda x: x.index.tolist())
            .to_dict()
        )
        first_neighbours = clean_dict(first_neighbours)

        # Create merged administrative regions
        merged_shapes = gpd.GeoDataFrame(
            columns=admin_shapes.columns,
            index=sorted(first_neighbours.keys()),
            crs=admin_shapes.crs,
        )
        merged_shapes["contains"] = admin_shapes.loc[merged_shapes.index, "contains"]
        merged_shapes = merged_shapes.apply(
            lambda row: create_merged_admin_region(
                row,
                first_neighbours,
                admin_shapes,
            ),
            axis=1,
        )
        merged_shapes.drop_duplicates(subset="geometry", keep="first", inplace=True)

        # Remove merged regions from admin_shapes and update admin_shapes
        list_merged_adm = list(set(chain.from_iterable(merged_shapes["contains"])))
        list_merged_adm = [adm for adm in list_merged_adm if adm in admin_shapes.index]
        admin_shapes = admin_shapes.drop(list_merged_adm)
        admin_shapes = pd.concat([admin_shapes, merged_shapes])

    return admin_shapes


def build_admin_shapes(
    n: pypsa.Network,
    nuts3_shapes: str,
    clustering: str,
    admin_levels: dict[str, int],
    countries: list[str],
) -> gpd.GeoDataFrame:
    """
    Builds administrative shapes for bus regions based on NUTS3 regions and custom clustering configuration.

    Parameters
    ----------
    nuts3_shapes (str) : Path to the file containing NUTS3 shapes.
    clustering (str) : Clustering method to use ('administrative' or other).
    admin_levels (dict[str, int]) : Dictionary mapping country codes to administrative levels.
    countries (list[str]) : List of country codes to include.

    Returns
    -------
    admin_shapes (gpd.GeoDataFrame): GeoDataFrame containing the administrative regions.
    """
    level_map = {
        0: "country",
        1: "level1",
        2: "level2",
        3: "level3",
        "bz": "bidding_zone",
    }

    adm1_countries = ["BA", "MD", "UA", "XK"]
    adm1_countries = list(set(adm1_countries).intersection(n.buses.country.unique()))

    nuts3_regions = gpd.read_file(nuts3_shapes)
    nuts3_regions = nuts3_regions.set_index(nuts3_regions.columns[0])

    level = admin_levels.get("level", 0)

    if clustering == "administrative":
        logger.info(f"Building bus regions at administrative level {level}")

        nuts3_regions["column"] = level_map[level]

        # Only keep the values whose keys are in countries
        country_level = {
            k: v for k, v in admin_levels.items() if (k != "level") and (k in countries)
        }
        if country_level:
            country_level_list = "\n".join(
                [f"- {k}: level {v}" for k, v in country_level.items()]
            )
            logger.info(
                f"Setting individual administrative levels for:\n{country_level_list}"
            )
            nuts3_regions.loc[
                nuts3_regions["country"].isin(country_level.keys()), "column"
            ] = (
                nuts3_regions.loc[
                    nuts3_regions["country"].isin(country_level.keys()), "country"
                ]
                .map(country_level)
                .map(level_map)
            )

        # If GB is in the countries, set the level, aggregate London area to level 1 due to converging issues
        if "GB" in countries and level != "bz":
            nuts3_regions.loc[nuts3_regions.level1 == "GBI", "column"] = "level1"

        nuts3_regions["admin"] = nuts3_regions.apply(
            lambda row: row[row["column"]], axis=1
        )
    else:
        logger.info("Building bus regions per country level.")
        nuts3_regions["admin"] = nuts3_regions["country"]

    # Group by busmap
    admin_shapes = nuts3_regions[["admin", "geometry", "country"]]
    admin_shapes = admin_shapes.dissolve("admin")

    # Identify regions that do not contain buses and merge them with smallest
    # neighbouring area of the same parent region (NUTS3 -> NUTS2 -> NUTS1 -> country)
    buses = n.buses[["x", "y", "country"]].copy()
    buses["geometry"] = gpd.points_from_xy(buses["x"], buses["y"])
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs=n.crs)

    # Obtain parents
    admin_shapes["parent"] = admin_shapes.index.str.replace("-", "").map(
        lambda x: x[:2] if len(x) == 2 else x[:-1]
    )
    admin_shapes.loc[admin_shapes.country.isin(adm1_countries), "parent"] = (
        admin_shapes.loc[admin_shapes.country.isin(adm1_countries), "country"]
    )

    # Check how many buses are in each region
    admin_shapes = gpd.sjoin(
        admin_shapes, buses[["country", "geometry"]], how="left", predicate="contains"
    )
    admin_shapes["substations"] = admin_shapes.apply(
        lambda x: 1 if x["country_left"] == x["country_right"] else 0, axis=1
    )
    admin_shapes.rename(columns={"country_left": "country"}, inplace=True)
    admin_shapes = admin_shapes.groupby(admin_shapes.index).agg(
        {
            "geometry": "first",
            "country": "first",
            "parent": "first",
            "substations": "sum",
        }
    )
    admin_shapes = gpd.GeoDataFrame(
        admin_shapes, geometry="geometry", crs=nuts3_regions.crs
    )
    admin_shapes["isempty"] = admin_shapes["substations"] == 0

    # Initiate contains column
    if "contains" not in admin_shapes.columns:
        admin_shapes["contains"] = admin_shapes.apply(lambda row: [row.name], axis=1)

    if admin_shapes["isempty"].any():
        logger.info(
            "Administrative regions without buses found. Merging with neighbouring regions of the same parent region."
        )

        admin_shapes = merge_regions_recursive(admin_shapes)

        # Find closest regions for remaining islands
        logger.info("Finding closest administrative regions for remaining islands.")
        b_island = admin_shapes["isempty"] == True
        admin_shapes.loc[b_island, "neighbours"] = admin_shapes.loc[b_island].apply(
            lambda row: [get_nearest_neighbour(row, admin_shapes)],
            axis=1,
        )
        # remove detached regions that don't contain substations
        admin_shapes = admin_shapes[~admin_shapes.isempty]
        admin_shapes = merge_regions_recursive(admin_shapes, neighbours_missing=False)

        # Update names
        admin_shapes.index.name = "admin"
        admin_shapes.reset_index(inplace=True)
        admin_shapes["admin"] = admin_shapes.apply(
            lambda row: update_names(row["contains"]),
            axis=1,
        )
        admin_shapes.set_index("admin", inplace=True)

    return admin_shapes[["country", "parent", "contains", "substations", "geometry"]]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("base_network")
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    mp.set_start_method("spawn", force=True)

    countries = snakemake.params.countries
    nuts3_shapes = snakemake.input.nuts3_shapes
    clustering = snakemake.params.get("clustering", "busmap")
    admin_levels = snakemake.params.get("admin_levels")

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

    admin_shapes = build_admin_shapes(
        n,
        nuts3_shapes,
        clustering,
        admin_levels,
        countries,
    )

    onshore_regions, offshore_regions, onshore_shapes, offshore_shapes = (
        build_bus_shapes(
            n,
            admin_shapes,
            offshore_shapes,
            countries,
        )
    )

    # Export network
    n.meta = snakemake.config
    n.export_to_netcdf(snakemake.output.base_network)

    # Export shapes
    onshore_shapes.to_file(snakemake.output.regions_onshore)
    # append_bus_shapes(n, shapes, "onshore")

    offshore_shapes.to_file(snakemake.output.regions_offshore)
    # append_bus_shapes(n, offshore_shapes, "offshore")

    # Convert contains columns into strings (pyogrio-friendly)
    admin_shapes["contains"] = admin_shapes["contains"].apply(lambda x: ",".join(x))
    admin_shapes["contains"] = admin_shapes["contains"].astype(str)
    admin_shapes.to_file(snakemake.output.admin_shapes)
    # append_bus_shapes(n, admin_shapes, "admin")
