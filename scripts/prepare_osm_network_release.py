# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import folium
import geopandas as gpd
import numpy as np
import pypsa
from _helpers import configure_logging, set_scenario_config
from base_network import _get_linetype_by_voltage
from shapely.wkt import loads

logger = logging.getLogger(__name__)


GEO_CRS = "EPSG:4326"
BUSES_COLUMNS = [
    "bus_id",
    "voltage",
    "dc",
    "symbol",
    "under_construction",
    "tags",
    "x",
    "y",
    "country",
    "geometry",
]
LINES_COLUMNS = [
    "line_id",
    "bus0",
    "bus1",
    "voltage",
    "i_nom",
    "circuits",
    "s_nom",
    "r",
    "x",
    "b",
    "length",
    "underground",
    "under_construction",
    "type",
    "tags",
    "geometry",
]
LINKS_COLUMNS = [
    "link_id",
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
TRANSFORMERS_COLUMNS = [
    "transformer_id",
    "bus0",
    "bus1",
    "voltage_bus0",
    "voltage_bus1",
    "s_nom",
    "geometry",
]
CONVERTERS_COLUMNS = [
    "converter_id",
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "geometry",
]


def export_clean_csv(df, columns, output_file):
    """
    Export a cleaned DataFrame to a CSV file.

    Args:
        df (pandas.DataFrame): The DataFrame to be exported.
        columns (list): A list of column names to include in the exported CSV file.
        output_file (str): The path to the output CSV file.

    Returns:
        None
    """
    rename_dict = {
        "Bus": "bus_id",
        "Line": "line_id",
        "Link": "link_id",
        "Transformer": "transformer_id",
        "v_nom": "voltage",
        "num_parallel": "circuits",
    }

    if "converter_id" in columns:
        rename_dict["Link"] = "converter_id"

    df.reset_index().rename(columns=rename_dict).loc[:, columns].replace(
        {True: "t", False: "f"}
    ).to_csv(output_file, index=False, quotechar="'")

    return None


def create_geometries(network, crs=GEO_CRS):
    """
    Create GeoDataFrames for different network components with specified coordinate reference system (CRS).

    Parameters
    ----------
        network (PyPSA Network): The network object containing buses, lines, links, converters, and transformers data.
        crs (str, optional): Coordinate reference system to be used for the GeoDataFrames. Defaults to GEO_CRS.

    Returns
    -------
    tuple: A tuple containing the following GeoDataFrames:
        - buses (GeoDataFrame): GeoDataFrame containing bus data with geometries.
        - lines (GeoDataFrame): GeoDataFrame containing line data with geometries.
        - links (GeoDataFrame): GeoDataFrame containing link data with geometries.
        - converters (GeoDataFrame): GeoDataFrame containing converter data with geometries.
        - transformers (GeoDataFrame): GeoDataFrame containing transformer data with geometries.
    """
    buses = network.buses.reset_index()[
        [
            "Bus",
            "v_nom",
            "dc",
            "symbol",
            "under_construction",
            "tags",
            "geometry",
        ]
    ]
    buses["geometry"] = buses.geometry.apply(lambda x: loads(x))
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs=crs)

    lines = network.lines.reset_index()[
        [
            "Line",
            "bus0",
            "bus1",
            "v_nom",
            "i_nom",
            "num_parallel",
            "s_nom",
            "r",
            "x",
            "b",
            "length",
            "underground",
            "under_construction",
            "type",
            "tags",
            "geometry",
        ]
    ]
    # Create shapely linestring from geometry column
    lines["geometry"] = lines.geometry.apply(lambda x: loads(x))
    lines = gpd.GeoDataFrame(lines, geometry="geometry", crs=crs)

    links = (
        network.links[~is_converter]
        .reset_index()
        .rename(columns={"voltage": "v_nom"})[
            [
                "Link",
                "bus0",
                "bus1",
                "v_nom",
                "p_nom",
                "length",
                "underground",
                "under_construction",
                "tags",
                "geometry",
            ]
        ]
    )
    links["geometry"] = links.geometry.apply(lambda x: loads(x))
    links = gpd.GeoDataFrame(links, geometry="geometry", crs=crs)

    converters = (
        network.links[is_converter]
        .reset_index()
        .rename(columns={"voltage": "v_nom"})[
            [
                "Link",
                "bus0",
                "bus1",
                "v_nom",
                "p_nom",
                "geometry",
            ]
        ]
    )
    converters["geometry"] = converters.geometry.apply(lambda x: loads(x))
    converters = gpd.GeoDataFrame(converters, geometry="geometry", crs=crs)

    transformers = network.transformers.reset_index()[
        [
            "Transformer",
            "bus0",
            "bus1",
            "voltage_bus0",
            "voltage_bus1",
            "s_nom",
            "geometry",
        ]
    ]
    transformers["geometry"] = transformers.geometry.apply(lambda x: loads(x))
    transformers = gpd.GeoDataFrame(transformers, geometry="geometry", crs=crs)

    return buses, lines, links, converters, transformers


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_osm_network_release")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Params
    line_types = snakemake.params.line_types

    network = pypsa.Network(snakemake.input.base_network)

    logger.info("Re-adding line types to network.")
    # Re-add line types
    network.lines.loc[:, "type"] = network.lines.v_nom.apply(
        lambda x: _get_linetype_by_voltage(x, line_types)
    )

    # Calculate dependent variables (r, x)
    logger.info("Calculating dependent variables for network (r, x).")
    network.calculate_dependent_values()

    # i_nom
    network.lines["i_nom"] = (
        (network.lines.s_nom / network.lines.v_nom / network.lines.num_parallel)
        .div(np.sqrt(3))
        .round(3)
    )  # kA

    # Rounding of dependent values
    network.lines.s_nom = network.lines.s_nom.round(3)
    network.lines.r = network.lines.r.round(6)
    network.lines.x = network.lines.x.round(6)
    network.lines.b = network.lines.b.round(8)

    # Convert v_nom and num_parallel to integers
    network.buses.v_nom = network.buses.v_nom.astype(int)
    network.lines.v_nom = network.lines.v_nom.astype(int)
    network.lines.num_parallel = network.lines.num_parallel.astype(int)
    network.links.voltage = network.links.voltage.astype(int)
    network.links.p_nom = network.links.p_nom.astype(int)
    network.transformers.voltage_bus0 = network.transformers.voltage_bus0.astype(int)
    network.transformers.voltage_bus1 = network.transformers.voltage_bus1.astype(int)
    network.transformers.s_nom = network.transformers.s_nom.astype(int)

    network.buses["dc"] = network.buses.pop("carrier").map({"DC": "t", "AC": "f"})
    network.lines.length = network.lines.length * 1e3
    network.links.length = network.links.length * 1e3

    # Sort alphabetically
    network.buses.sort_index(inplace=True)
    network.transformers.sort_index(inplace=True)
    network.lines.sort_index(inplace=True)
    network.links.sort_index(inplace=True)

    # Export to clean csv for release
    logger.info(f"Exporting {len(network.buses)} buses to %s", snakemake.output.buses)
    export_clean_csv(network.buses, BUSES_COLUMNS, snakemake.output.buses)

    logger.info(
        f"Exporting {len(network.transformers)} transformers to %s",
        snakemake.output.transformers,
    )
    export_clean_csv(
        network.transformers, TRANSFORMERS_COLUMNS, snakemake.output.transformers
    )

    logger.info(f"Exporting {len(network.lines)} lines to %s", snakemake.output.lines)
    export_clean_csv(network.lines, LINES_COLUMNS, snakemake.output.lines)

    # Boolean that specifies if link element is a converter
    is_converter = network.links.index.str.startswith("conv") == True

    logger.info(
        f"Exporting {len(network.links[~is_converter])} links to %s",
        snakemake.output.links,
    )
    export_clean_csv(
        network.links[~is_converter], LINKS_COLUMNS, snakemake.output.links
    )

    logger.info(
        f"Exporting {len(network.links[is_converter])} converters to %s",
        snakemake.output.converters,
    )
    export_clean_csv(
        network.links[is_converter], CONVERTERS_COLUMNS, snakemake.output.converters
    )

    ### Create interactive map
    buses, lines, links, converters, transformers = create_geometries(
        network, crs=GEO_CRS
    )
    stations_polygon = gpd.read_file(snakemake.input.stations_polygon)
    buses_polygon = gpd.read_file(snakemake.input.buses_polygon)

    # Only keep stations_polygon that contain buses points
    stations_polygon = gpd.sjoin(
        stations_polygon, buses, how="left", predicate="contains"
    )
    stations_polygon = stations_polygon[stations_polygon.index_right.notnull()]
    stations_polygon = stations_polygon.drop_duplicates(subset=["station_id"])
    stations_polygon = stations_polygon[["station_id", "geometry"]]

    buses_polygon = gpd.sjoin(buses_polygon, buses, how="left", predicate="contains")
    buses_polygon = buses_polygon[buses_polygon.index_right.notnull()]
    buses_polygon = buses_polygon.drop_duplicates(subset=["bus_id", "dc_left"])
    buses_polygon.rename(columns={"dc_left": "dc"}, inplace=True)
    buses_polygon = buses_polygon[["bus_id", "dc", "geometry"]]

    map = None
    map = folium.Map(tiles="CartoDB positron", zoom_start=5, location=[53.5, 10])
    map = stations_polygon.loc[
        (
            stations_polygon.station_id.str.startswith("way")
            | stations_polygon.station_id.str.startswith("relation")
        )
    ].explore(
        color="darkred",
        popup=True,
        m=map,
        name="Clustered substations",
        zindex=100,
    )
    map = stations_polygon.loc[
        ~(
            stations_polygon.station_id.str.startswith("way")
            | stations_polygon.station_id.str.startswith("relation")
        )
    ].explore(
        color="grey",
        popup=True,
        m=map,
        name="Clustered substations (virtual)",
        zindex=101,
    )
    map = buses_polygon.loc[buses_polygon.dc == False].explore(
        color="yellow",
        popup=True,
        m=map,
        name="Buses (AC)",
        zindex=102,
    )
    map = buses_polygon.loc[buses_polygon.dc == True].explore(
        color="grey",
        popup=True,
        m=map,
        name="Buses (DC)",
        zindex=103,
    )
    map = lines.explore(
        color="rosybrown",
        popup=True,
        m=map,
        name="Lines (AC)",
        zindex=104,
    )
    map = links.explore(
        color="darkseagreen",
        popup=True,
        m=map,
        name="Links (DC)",
        zindex=105,
    )
    map = transformers.explore(
        color="orange",
        popup=True,
        m=map,
        name="Transformers",
        zindex=106,
    )
    map = converters.explore(
        color="purple",
        popup=True,
        m=map,
        name="Converters",
        zindex=107,
    )
    map = buses.loc[buses.dc == "f"].explore(
        color="red",
        popup=True,
        m=map,
        name="Buses (AC, Points)",
        zindex=108,
    )
    map = buses.loc[buses.dc == "t"].explore(
        color="black",
        popup=True,
        m=map,
        name="Buses (DC, Points)",
        zindex=109,
    )
    # Add legend
    folium.LayerControl(collapsed=False).add_to(map)

    map_title = "Prebuilt electricity high-voltage grid based on OpenStreetMap data"
    map.get_root().html.add_child(
        folium.Element(
            f"<h4 style='position:absolute;z-index:100000;left:1vw;bottom:5px'>{map_title}</h4>"
        )
    )
    map

    # Export map
    logger.info("Exporting interactive map.")
    map.save(snakemake.output.map)

    logger.info("Export of OSM network for release complete.")
