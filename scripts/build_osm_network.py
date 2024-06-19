# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur and PyPSA-Earth Authors
#
# SPDX-License-Identifier: MIT

import logging
import os

import geopandas as gpd
import numpy as np
import pandas as pd
from _benchmark import memory_logger
from _helpers import configure_logging, set_scenario_config
from shapely.geometry import LineString, Point
from shapely.ops import linemerge, split
from tqdm import tqdm

logger = logging.getLogger(__name__)

# list of recognised nan values (NA and na excluded as may be confused with Namibia 2-letter country code)
NA_VALUES = ["NULL", "", "N/A", "NAN", "NaN", "nan", "Nan", "n/a", "null"]


def read_csv_nafix(file, **kwargs):
    "Function to open a csv as pandas file and standardize the na value"
    if "keep_default_na" not in kwargs:
        kwargs["keep_default_na"] = False
    if "na_values" not in kwargs:
        kwargs["na_values"] = NA_VALUES

    if os.stat(file).st_size > 0:
        return pd.read_csv(file, **kwargs)
    else:
        return pd.DataFrame()


def save_to_geojson(df, fn):
    """
    Save a (Geo)DataFrame to a GeoJSON file.

    Parameters:
    - df: The (Geo)DataFrame to be saved.
    - fn: The filename (including the path) of the output GeoJSON file.

    Returns:
    None
    """
    if os.path.exists(fn):
        os.unlink(fn)  # remove file if it exists

    # save file if the (Geo)DataFrame is non-empty
    if df.empty:
        # create empty file to avoid issues with snakemake
        with open(fn, "w") as fp:
            pass
    else:
        # save file
        df.to_file(fn, driver="GeoJSON")


def read_geojson(fn, cols=[], dtype=None, crs="EPSG:4326"):
    """
    Function to read a geojson file fn. When the file is empty, then an empty
    GeoDataFrame is returned having columns cols, the specified crs and the
    columns specified by the dtype dictionary it not none.

    Parameters:
    ------------
    fn : str
        Path to the file to read
    cols : list
        List of columns of the GeoDataFrame
    dtype : dict
        Dictionary of the type of the object by column
    crs : str
        CRS of the GeoDataFrame
    """
    # if the file is non-zero, read the geodataframe and return it
    if os.path.getsize(fn) > 0:
        return gpd.read_file(fn)
    else:
        # else return an empty GeoDataFrame
        df = gpd.GeoDataFrame(columns=cols, geometry=[], crs=crs)
        if isinstance(dtype, dict):
            for k, v in dtype.items():
                df[k] = df[k].astype(v)
        return df


def to_csv_nafix(df, path, **kwargs):
    """
    Write a pandas DataFrame to a CSV file with NA values replaced.

    Parameters:
    - df: pandas DataFrame
        The DataFrame to be written to the CSV file.
    - path: str
        The file path where the CSV file will be saved.
    - **kwargs: keyword arguments
        Additional arguments to be passed to the `to_csv` function of pandas.

    Returns:
    - None

    If the DataFrame is not empty or does not have empty columns, it will be
    written to the CSV file with NA values replaced by the first value in the
    `NA_VALUES` list. If the DataFrame is empty or has empty columns, an empty
    file will be created at the specified path.
    """
    if "na_rep" in kwargs:
        del kwargs["na_rep"]
    if not df.empty or not df.columns.empty:
        return df.to_csv(path, **kwargs, na_rep=NA_VALUES[0])
    else:
        with open(path, "w") as fp:
            pass


def line_endings_to_bus_conversion(lines):
    """
    Converts line endings to bus connections.

    This function takes a df of lines and converts the line endings to bus
    connections. It performs the necessary operations to ensure that the line
    endings are properly connected to the buses in the network.

    Parameters:
    lines (DataFrame)

    Returns:
    lines (DataFrame)
    """
    # Assign to every line a start and end point

    lines["bounds"] = lines["geometry"].boundary  # create start and end point

    lines["bus_0_coors"] = lines["bounds"].map(lambda p: p.geoms[0])
    lines["bus_1_coors"] = lines["bounds"].map(lambda p: p.geoms[-1])

    # splits into coordinates
    lines["bus0_lon"] = lines["bus_0_coors"].x
    lines["bus0_lat"] = lines["bus_0_coors"].y
    lines["bus1_lon"] = lines["bus_1_coors"].x
    lines["bus1_lat"] = lines["bus_1_coors"].y

    return lines


# tol in m
def set_substations_ids(buses, distance_crs, tol=5000):
    """
    Function to set substations ids to buses, accounting for location
    tolerance.

    The algorithm is as follows:

    1. initialize all substation ids to -1
    2. if the current substation has been already visited [substation_id < 0], then skip the calculation
    3. otherwise:
        1. identify the substations within the specified tolerance (tol)
        2. when all the substations in tolerance have substation_id < 0, then specify a new substation_id
        3. otherwise, if one of the substation in tolerance has a substation_id >= 0, then set that substation_id to all the others;
           in case of multiple substations with substation_ids >= 0, the first value is picked for all
    """

    buses["station_id"] = -1

    # create temporary series to execute distance calculations using m as reference distances
    temp_bus_geom = buses.geometry.to_crs(distance_crs)

    # set tqdm options for substation ids
    tqdm_kwargs_substation_ids = dict(
        ascii=False,
        unit=" buses",
        total=buses.shape[0],
        desc="Set substation ids ",
    )

    station_id = 0
    for i, row in tqdm(buses.iterrows(), **tqdm_kwargs_substation_ids):
        if buses.loc[i, "station_id"] >= 0:
            continue

        # get substations within tolerance
        close_nodes = np.flatnonzero(
            temp_bus_geom.distance(temp_bus_geom.loc[i]) <= tol
        )

        if len(close_nodes) == 1:
            # if only one substation is in tolerance, then the substation is the current one iì
            # Note that the node cannot be with substation_id >= 0, given the preliminary check
            # at the beginning of the for loop
            buses.loc[buses.index[i], "station_id"] = station_id
            # update station id
            station_id += 1
        else:
            # several substations in tolerance
            # get their ids
            subset_substation_ids = buses.loc[buses.index[close_nodes], "station_id"]
            # check if all substation_ids are negative (<0)
            all_neg = subset_substation_ids.max() < 0
            # check if at least a substation_id is negative (<0)
            some_neg = subset_substation_ids.min() < 0

            if all_neg:
                # when all substation_ids are negative, then this is a new substation id
                # set the current station_id and increment the counter
                buses.loc[buses.index[close_nodes], "station_id"] = station_id
                station_id += 1
            elif some_neg:
                # otherwise, when at least a substation_id is non-negative, then pick the first value
                # and set it to all the other substations within tolerance
                sub_id = -1
                for substation_id in subset_substation_ids:
                    if substation_id >= 0:
                        sub_id = substation_id
                        break
                buses.loc[buses.index[close_nodes], "station_id"] = sub_id


def set_lines_ids(lines, buses, distance_crs):
    """
    Function to set line buses ids to the closest bus in the list.
    """
    # set tqdm options for set lines ids
    tqdm_kwargs_line_ids = dict(
        ascii=False,
        unit=" lines",
        total=lines.shape[0],
        desc="Set line bus ids ",
    )

    # initialization
    lines["bus0"] = -1
    lines["bus1"] = -1

    busesepsg = buses.to_crs(distance_crs)
    linesepsg = lines.to_crs(distance_crs)

    for i, row in tqdm(linesepsg.iterrows(), **tqdm_kwargs_line_ids):
        # select buses having the voltage level of the current line
        buses_sel = busesepsg[
            (buses["voltage"] == row["voltage"]) & (buses["dc"] == row["dc"])
        ]

        # find the closest node of the bus0 of the line
        bus0_id = buses_sel.geometry.distance(row.geometry.boundary.geoms[0]).idxmin()
        lines.loc[i, "bus0"] = buses.loc[bus0_id, "bus_id"]

        # check if the line starts exactly in the node, otherwise modify the linestring
        distance_bus0 = busesepsg.geometry.loc[bus0_id].distance(
            row.geometry.boundary.geoms[0]
        )
        if distance_bus0 > 0.0:
            # the line does not start in the node, thus modify the linestring
            lines.loc[i, "geometry"] = linemerge(
                [
                    LineString(
                        [
                            buses.geometry.loc[bus0_id],
                            lines.geometry.loc[i].boundary.geoms[0],
                        ]
                    ),
                    lines.geometry.loc[i],
                ]
            )

        # find the closest node of the bus1 of the line
        bus1_id = buses_sel.geometry.distance(row.geometry.boundary.geoms[1]).idxmin()
        lines.loc[i, "bus1"] = buses.loc[bus1_id, "bus_id"]

        # check if the line ends exactly in the node, otherwise modify the linestring
        distance_bus1 = busesepsg.geometry.loc[bus1_id].distance(
            row.geometry.boundary.geoms[1]
        )
        if distance_bus1 > 0.0:
            # the line does not end in the node, thus modify the linestring
            lines.loc[i, "geometry"] = linemerge(
                [
                    lines.geometry.loc[i],
                    LineString(
                        [
                            lines.geometry.loc[i].boundary.geoms[1],
                            buses.geometry.loc[bus1_id],
                        ]
                    ),
                ]
            )

    return lines, buses


def merge_stations_same_station_id(
    buses, delta_lon=0.001, delta_lat=0.001, precision=4
):
    """
    Function to merge buses with same voltage and station_id This function
    iterates over all substation ids and creates a bus_id for every substation
    and voltage level.

    Therefore, a substation with multiple voltage levels is represented
    with different buses, one per voltage level
    """
    # initialize list of cleaned buses
    buses_clean = []

    # initialize the number of buses
    n_buses = 0

    for g_name, g_value in buses.groupby(by="station_id"):
        # average location of the buses having the same station_id
        station_point_x = np.round(g_value.geometry.x.mean(), precision)
        station_point_y = np.round(g_value.geometry.y.mean(), precision)
        # is_dclink_boundary_point = any(g_value["is_dclink_boundary_point"])

        # loop for every voltage level in the bus
        # The location of the buses is averaged; in the case of multiple voltage levels for the same station_id,
        # each bus corresponding to a voltage level and each polatity is located at a distance regulated by delta_lon/delta_lat
        v_it = 0
        for v_name, bus_row in g_value.groupby(by=["voltage", "dc"]):
            lon_bus = np.round(station_point_x + v_it * delta_lon, precision)
            lat_bus = np.round(station_point_y + v_it * delta_lat, precision)

            # add the bus
            buses_clean.append(
                [
                    n_buses,  # "bus_id"
                    g_name,  # "station_id"
                    v_name[0],  # "voltage"
                    bus_row["dc"].all(),  # "dc"
                    "|".join(bus_row["symbol"].unique()),  # "symbol"
                    bus_row["under_construction"].any(),  # "under_construction"
                    "|".join(bus_row["tag_substation"].unique()),  # "tag_substation"
                    bus_row["tag_area"].sum(),  # "tag_area"
                    lon_bus,  # "lon"
                    lat_bus,  # "lat"
                    bus_row["country"].iloc[0],  # "country",
                    # is_dclink_boundary_point,  # check if new bus was formed of at least one DC link boundary point
                    Point(
                        lon_bus,
                        lat_bus,
                    ),  # "geometry"
                ]
            )

            # increase counters
            v_it += 1
            n_buses += 1

    # names of the columns
    buses_clean_columns = [
        "bus_id",
        "station_id",
        "voltage",
        "dc",
        "symbol",
        "under_construction",
        "tag_substation",
        "tag_area",
        "x",
        "y",
        "country",
        # "is_dclink_boundary_point",
        "geometry",
    ]

    gdf_buses_clean = gpd.GeoDataFrame(
        buses_clean, columns=buses_clean_columns
    ).set_crs(crs=buses.crs, inplace=True)

    return gdf_buses_clean


def get_ac_frequency(df, fr_col="tag_frequency"):
    """
    # Function to define a default frequency value.

    Attempts to find the most usual non-zero frequency across the
    dataframe; 50 Hz is assumed as a back-up value
    """

    # Initialize a default frequency value
    ac_freq_default = 50

    grid_freq_levels = df[fr_col].value_counts(sort=True, dropna=True)
    if not grid_freq_levels.empty:
        # AC lines frequency shouldn't be 0Hz
        ac_freq_levels = grid_freq_levels.loc[
            grid_freq_levels.index.get_level_values(0) != "0"
        ]
        ac_freq_default = ac_freq_levels.index.get_level_values(0)[0]

    return ac_freq_default


def get_transformers(buses, lines):
    """
    Function to create fake transformer lines that connect buses of the same
    station_id at different voltage.
    """

    ac_freq = get_ac_frequency(lines)
    df_transformers = []

    # Transformers should be added between AC buses only
    # TODO pypsa-eur: Fix this! instead of tilde use !=
    buses_ac = buses[buses["dc"] != True]
    for g_name, g_value in buses_ac.sort_values("voltage", ascending=True).groupby(
        by="station_id"
    ):
        # note: by construction there cannot be more that two buses with the same station_id and same voltage
        n_voltages = len(g_value)

        if n_voltages > 1:
            for id in range(0, n_voltages - 1):
                # when g_value has more than one node, it means that there are multiple voltages for the same bus
                geom_trans = LineString(
                    [g_value.geometry.iloc[id], g_value.geometry.iloc[id + 1]]
                )

                df_transformers.append(
                    [
                        f"transf_{g_name}_{id}",  # "line_id"
                        g_value["bus_id"].iloc[id],  # "bus0"
                        g_value["bus_id"].iloc[id + 1],  # "bus1"
                        g_value.voltage.iloc[id],  # "voltage_bus0"
                        g_value.voltage.iloc[id + 1],  # "voltage_bus0"
                        g_value.country.iloc[id],  # "country"
                        geom_trans,  # "geometry"
                    ]
                )
    # TODO pypsa-eur: fix bug in pypsa-earth, where the id column is wrongly named "line_id" instead of "transformer_id
    # name of the columns
    trasf_columns = [
        "transformer_id",
        "bus0",
        "bus1",
        "voltage_bus0",
        "voltage_bus1",
        "country",
        "geometry",
    ]

    df_transformers = gpd.GeoDataFrame(df_transformers, columns=trasf_columns)
    if not df_transformers.empty:
        init_index = 0 if lines.empty else lines.index[-1] + 1
        df_transformers.set_index(init_index + df_transformers.index, inplace=True)
    # update line endings
    df_transformers = line_endings_to_bus_conversion(df_transformers)

    return df_transformers


def get_converters(buses):
    """
    Function to create fake converter lines that connect buses of the same
    station_id of different polarities.
    """

    df_converters = []

    for g_name, g_value in buses.sort_values("voltage", ascending=True).groupby(
        by="station_id"
    ):
        # note: by construction there cannot be more that two buses with the same station_id and same voltage
        n_voltages = len(g_value)

        # A converter stations should have both AC and DC parts
        if g_value["dc"].any() & ~g_value["dc"].all():
            dc_voltage = g_value[g_value.dc]["voltage"].values

            for u in dc_voltage:
                id_0 = g_value[g_value["dc"] & g_value["voltage"].isin([u])].index[0]

                ac_voltages = g_value[~g_value.dc]["voltage"]
                # A converter is added between a DC nodes and AC one with the closest voltage
                id_1 = ac_voltages.sub(u).abs().idxmin()

                geom_conv = LineString(
                    [g_value.geometry.loc[id_0], g_value.geometry.loc[id_1]]
                )

                # check if bus is a dclink boundary point, only then add converter
                df_converters.append(
                    [
                        f"convert_{g_name}_{id_0}",  # "line_id"
                        g_value["bus_id"].loc[id_0],  # "bus0"
                        g_value["bus_id"].loc[id_1],  # "bus1"
                        False,  # "underground"
                        False,  # "under_construction"
                        g_value.country.loc[id_0],  # "country"
                        geom_conv,  # "geometry"
                    ]
                )

    # name of the columns
    conv_columns = [
        "converter_id",
        "bus0",
        "bus1",
        "underground",
        "under_construction",
        "country",
        "geometry",
    ]

    df_converters = gpd.GeoDataFrame(df_converters, columns=conv_columns).reset_index()

    return df_converters


def connect_stations_same_station_id(lines, buses):
    """
    Function to create fake links between substations with the same
    substation_id.
    """
    ac_freq = get_ac_frequency(lines)
    station_id_list = buses.station_id.unique()

    add_lines = []
    from shapely.geometry import LineString

    for s_id in station_id_list:
        buses_station_id = buses[buses.station_id == s_id]

        if len(buses_station_id) > 1:
            for b_it in range(1, len(buses_station_id)):
                add_lines.append(
                    [
                        f"link{buses_station_id}_{b_it}",  # "line_id"
                        buses_station_id.index[0],  # "bus0"
                        buses_station_id.index[b_it],  # "bus1"
                        400000,  # "voltage"
                        1,  # "circuits"
                        0.0,  # "length"
                        False,  # "underground"
                        False,  # "under_construction"
                        "transmission",  # "tag_type"
                        ac_freq,  # "tag_frequency"
                        buses_station_id.country.iloc[0],  # "country"
                        LineString(
                            [
                                buses_station_id.geometry.iloc[0],
                                buses_station_id.geometry.iloc[b_it],
                            ]
                        ),  # "geometry"
                        LineString(
                            [
                                buses_station_id.geometry.iloc[0],
                                buses_station_id.geometry.iloc[b_it],
                            ]
                        ).bounds,  # "bounds"
                        buses_station_id.geometry.iloc[0],  # "bus_0_coors"
                        buses_station_id.geometry.iloc[b_it],  # "bus_1_coors"
                        buses_station_id.lon.iloc[0],  # "bus0_lon"
                        buses_station_id.lat.iloc[0],  # "bus0_lat"
                        buses_station_id.lon.iloc[b_it],  # "bus1_lon"
                        buses_station_id.lat.iloc[b_it],  # "bus1_lat"
                    ]
                )

    # name of the columns
    add_lines_columns = [
        "line_id",
        "bus0",
        "bus1",
        "voltage",
        "circuits",
        "length",
        "underground",
        "under_construction",
        "tag_type",
        "tag_frequency",
        "country",
        "geometry",
        "bounds",
        "bus_0_coors",
        "bus_1_coors",
        "bus0_lon",
        "bus0_lat",
        "bus1_lon",
        "bus1_lat",
    ]

    df_add_lines = gpd.GeoDataFrame(pd.concat(add_lines), columns=add_lines_columns)
    lines = pd.concat([lines, df_add_lines], ignore_index=True)

    return lines


def set_lv_substations(buses):
    """
    Function to set what nodes are lv, thereby setting substation_lv The
    current methodology is to set lv nodes to buses where multiple voltage
    level are found, hence when the station_id is duplicated.
    """
    # initialize column substation_lv to true
    buses["substation_lv"] = True

    # For each station number with multiple buses make lowest voltage `substation_lv = TRUE`
    bus_with_stations_duplicates = buses[
        buses.station_id.duplicated(keep=False)
    ].sort_values(by=["station_id", "voltage"])
    lv_bus_at_station_duplicates = (
        buses[buses.station_id.duplicated(keep=False)]
        .sort_values(by=["station_id", "voltage"])
        .drop_duplicates(subset=["station_id"])
    )
    # Set all buses with station duplicates "False"
    buses.loc[bus_with_stations_duplicates.index, "substation_lv"] = False
    # Set lv_buses with station duplicates "True"
    buses.loc[lv_bus_at_station_duplicates.index, "substation_lv"] = True

    return buses


def merge_stations_lines_by_station_id_and_voltage(
    lines, links, buses, distance_crs, tol=5000
):
    """
    Function to merge close stations and adapt the line datasets to adhere to
    the merged dataset.
    """

    logger.info(" - Setting substation ids with tolerance of %.2f m" % (tol))

    # TODO pypsa-eur: Add this fix to pypsa-earth: Buses should not be clustered geographically if they are different
    # bus types (AC != DC)
    buses_ac = buses[buses["dc"] == False].reset_index()
    buses_dc = buses[buses["dc"] == True].reset_index()

    # set substation ids
    # set_substations_ids(buses, distance_crs, tol=tol)
    set_substations_ids(buses_ac, distance_crs, tol=tol)
    set_substations_ids(buses_dc, distance_crs, tol=tol)

    logger.info(" - Merging substations with the same id")

    # merge buses with same station id and voltage
    if not buses.empty:
        buses_ac = merge_stations_same_station_id(buses_ac)
        buses_dc = merge_stations_same_station_id(buses_dc)
        buses_dc["bus_id"] = buses_ac["bus_id"].max() + buses_dc["bus_id"] + 1
        buses = pd.concat([buses_ac, buses_dc], ignore_index=True)
        set_substations_ids(buses, distance_crs, tol=tol)

    logger.info(" - Specifying the bus ids of the line endings")

    # set the bus ids to the line dataset
    lines, buses = set_lines_ids(lines, buses, distance_crs)
    links, buses = set_lines_ids(links, buses, distance_crs)

    # drop lines starting and ending in the same node
    lines.drop(lines[lines["bus0"] == lines["bus1"]].index, inplace=True)
    links.drop(links[links["bus0"] == links["bus1"]].index, inplace=True)
    # update line endings
    lines = line_endings_to_bus_conversion(lines)
    links = line_endings_to_bus_conversion(links)

    # set substation_lv
    set_lv_substations(buses)

    # reset index
    lines.reset_index(drop=True, inplace=True)
    links.reset_index(drop=True, inplace=True)
    # if len(links) > 0:
    #     links.reset_index(drop=True, inplace=True)

    return lines, links, buses


def build_network(
    inputs,
    outputs,
    geo_crs,
    distance_crs,
):
    osm_clean_columns = {
        "substation": {
            "bus_id": "object",
            "station_id": "float",
            "voltage": "float",
            "dc": "bool",
            "symbol": "object",
            "under_construction": "bool",
            "tag_substation": "str",
            "tag_area": "str",
            "lon": "float",
            "lat": "float",
            "country": "str",
            "geometry": "object",
            "tag_source": "str",
        },
        "line": {
            "line_id": "object",
            "bus0": "object",
            "bus1": "object",
            "voltage": "float",
            "circuits": "float",
            "length": "float",
            "underground": "bool",
            "under_construction": "bool",
            "tag_type": "str",
            "tag_frequency": "float",
            "dc": "bool",
            "country": "object",
            "geometry": "object",
        },
        "link": {
            "link_id": "object",
            "bus0": "object",
            "bus1": "object",
            "voltage": "float",
            "length": "float",
            "under_construction": "bool",
            "dc": "bool",
            "country": "object",
            "geometry": "object",
        },
    }

    logger.info("Reading input data.")
    buses = read_geojson(
        inputs["substations"],
        osm_clean_columns["substation"].keys(),
        dtype=osm_clean_columns["substation"],
    )

    lines = read_geojson(
        inputs["lines"],
        osm_clean_columns["line"].keys(),
        dtype=osm_clean_columns["line"],
    )

    links = read_geojson(
        inputs["links"],
        osm_clean_columns["link"].keys(),
        dtype=osm_clean_columns["link"],
    )

    lines = line_endings_to_bus_conversion(lines)
    links = line_endings_to_bus_conversion(links)

    # METHOD to merge buses with same voltage and within tolerance
    tol = snakemake.config["electricity_network"]["osm_group_tolerance_buses"]
    logger.info(f"Aggregating close substations: Enabled with tolerance {tol} m")

    lines, links, buses = merge_stations_lines_by_station_id_and_voltage(
        lines, links, buses, distance_crs, tol=tol
    )

    # Recalculate lengths of lines
    utm = lines.estimate_utm_crs(datum_name="WGS 84")
    lines["length"] = lines.to_crs(utm).length
    links["length"] = links.to_crs(utm).length

    # TODO pypsa-eur: check if needed for updated links scripts
    # get transformers: modelled as lines connecting buses with different voltage
    transformers = get_transformers(buses, lines)

    # get converters: currently modelled as links connecting buses with different polarity
    converters = get_converters(buses)

    logger.info("Saving outputs")

    # create clean directory if not already exist
    if not os.path.exists(outputs["lines"]):
        os.makedirs(os.path.dirname(outputs["lines"]), exist_ok=True)

    ### Convert output to pypsa-eur friendly format
    # Rename "substation" in buses["symbol"] to "Substation"
    buses["symbol"] = buses["symbol"].replace({"substation": "Substation"})

    # Drop unnecessary index column and set respective element ids as index
    lines.set_index("line_id", inplace=True)
    links.set_index("link_id", inplace=True)
    converters.set_index("converter_id", inplace=True)
    transformers.set_index("transformer_id", inplace=True)
    buses.set_index("bus_id", inplace=True)

    # Convert voltages from V to kV
    lines["voltage"] = lines["voltage"] / 1000
    links["voltage"] = links["voltage"] / 1000
    transformers["voltage_bus0"], transformers["voltage_bus1"] = (
        transformers["voltage_bus0"] / 1000,
        transformers["voltage_bus1"] / 1000,
    )
    buses["voltage"] = buses["voltage"] / 1000

    # Convert 'true' and 'false' to 't' and 'f'
    lines = lines.replace({True: "t", False: "f"})
    links = links.replace({True: "t", False: "f"})
    converters = converters.replace({True: "t", False: "f"})
    buses = buses.replace({True: "t", False: "f"})

    # Change column orders
    cols_lines = [
        "bus0",
        "bus1",
        "voltage",
        "circuits",
        "length",
        "underground",
        "under_construction",
        "geometry",
    ]

    lines = lines[cols_lines]

    cols_links = [
        "bus0",
        "bus1",
        "voltage",
        "p_nom",
        "length",
        "under_construction",
        "geometry",
    ]

    links = links[cols_links]

    cols_transformers = [
        "bus0",
        "bus1",
        "voltage_bus0",
        "voltage_bus1",
        "country",
        "geometry",
    ]

    transformers = transformers[cols_transformers]

    to_csv_nafix(lines, outputs["lines"], quotechar="'")  # Generate CSV
    to_csv_nafix(links, outputs["links"], quotechar="'")  # Generate CSV
    to_csv_nafix(converters, outputs["converters"], quotechar="'")  # Generate CSV
    to_csv_nafix(transformers, outputs["transformers"], quotechar="'")  # Generate CSV

    # Export to GeoJSON for quick validations
    save_to_geojson(
        gpd.GeoDataFrame(lines),
        outputs["lines_geojson"],
    )
    save_to_geojson(
        gpd.GeoDataFrame(links),
        outputs["links_geojson"],
    )
    save_to_geojson(
        gpd.GeoDataFrame(converters, geometry="geometry", crs=geo_crs),
        outputs["converters_geojson"],
    )
    save_to_geojson(
        gpd.GeoDataFrame(transformers, geometry="geometry", crs=geo_crs),
        outputs["transformers_geojson"],
    )

    # create clean directory if not already exist
    if not os.path.exists(outputs["substations"]):
        os.makedirs(os.path.dirname(outputs["substations"]), exist_ok=True)
    # Generate CSV
    to_csv_nafix(buses, outputs["substations"], quotechar="'")
    save_to_geojson(
        gpd.GeoDataFrame(buses, geometry="geometry", crs=geo_crs),
        outputs["substations_geojson"],
    )

    return None


# Function to check if two lines are connected
def are_lines_connected(line1, line2):
    """
    Check if two lines are connected.

    Parameters:
    line1 (dict): A dictionary representing the first line.
    line2 (dict): A dictionary representing the second line.

    Returns:
    tuple: A tuple of boolean values indicating the connection status between
    the lines.

    The tuple contains four elements:
    - True if the first line's bus_0_coors is almost equal to the second line's
      bus_0_coors, False otherwise.
    - True if the first line's bus_0_coors is almost equal to the second line's
      bus_1_coors, False otherwise.
    - True if the first line's bus_1_coors is almost equal to the second line's
      bus_0_coors, False otherwise.
    - True if the first line's bus_1_coors is almost equal to the second line's
      bus_1_coors, False otherwise.
    """
    return (
        are_almost_equal(line1["bus_0_coors"], line2["bus_0_coors"]),
        are_almost_equal(line1["bus_0_coors"], line2["bus_1_coors"]),
        are_almost_equal(line1["bus_1_coors"], line2["bus_0_coors"]),
        are_almost_equal(line1["bus_1_coors"], line2["bus_1_coors"]),
    )


def _dfs(adj_matrix, visited, current_vertex, path):
    """
    Perform a depth-first search (DFS) on a graph represented by an adjacency
    matrix.

    Parameters:
    - adj_matrix (list of lists): The adjacency matrix representing the graph.
    - visited (list of bool): A list to keep track of visited vertices.
    - current_vertex (int): The current vertex being visited.
    - path (list): The path of vertices visited so far.

    Returns:
    - path (list): The path of vertices visited during the DFS.
    """
    visited[current_vertex] = True
    path.append(current_vertex)
    for neighbor in range(len(adj_matrix)):
        if adj_matrix[current_vertex][neighbor] == 1 and not visited[neighbor]:
            _dfs(adj_matrix, visited, neighbor, path)
    return path


# Returns all connected paths as a vector
def find_paths(adj_matrix):
    """
    Find all paths in a graph represented by an adjacency matrix.

    Parameters:
    - adj_matrix (list of lists): The adjacency matrix representing the graph.

    Returns:
    - paths (list of lists): A list of lists, where each inner list represents
      a path in the graph.
    """
    visited = [False] * len(adj_matrix)
    paths = []
    for vertex in range(len(adj_matrix)):
        if not visited[vertex]:
            path = _dfs(adj_matrix, visited, vertex, [])
            if path:
                paths.append(path)
    return paths


def are_almost_equal(point1, point2, tolerance=1e-6):
    """
    Check if two Shapely points are almost equal with a given tolerance.

    Args:
    point1 (Point): First Shapely point.
    point2 (Point): Second Shapely point.
    tolerance (float): Tolerance for coordinate deviation.

    Returns:
    bool: True if the points are almost equal, False otherwise.
    """
    return abs(point1.x - point2.x) < tolerance and abs(point1.y - point2.y) < tolerance


def merge_linestrings(gdf):
    """
    Merge LineStrings in a GeoDataFrame wherever the endpoints match.

    Parameters:
    gdf (GeoDataFrame): A GeoDataFrame containing LineString geometries.

    Returns:
    GeoDataFrame: A GeoDataFrame with merged LineString geometries.
    """
    gdf = gdf.copy()
    if len(gdf) == 1:
        return gdf

    lines = list(gdf.geometry)
    merged_lines = []
    while lines:
        line = lines.pop(0)
        merged_line = line
        i = 0
        while i < len(lines):
            if are_almost_equal(
                Point(merged_line.coords[-1]), Point(lines[i].coords[0])
            ):
                merged_line = LineString(
                    list(merged_line.coords) + list(lines.pop(i).coords[1:])
                )
                i = 0  # Restart the scan after merging
            elif are_almost_equal(
                Point(merged_line.coords[0]), Point(lines[i].coords[-1])
            ):
                merged_line = LineString(
                    list(lines.pop(i).coords)[:-1] + list(merged_line.coords)
                )
                i = 0  # Restart the scan after merging
            elif are_almost_equal(
                Point(merged_line.coords[-1]), Point(lines[i].coords[-1])
            ):
                merged_line = LineString(
                    list(merged_line.coords) + list(lines.pop(i).coords[::-1])[1:]
                )
                i = 0  # Restart the scan after merging
            elif are_almost_equal(
                Point(merged_line.coords[0]), Point(lines[i].coords[0])
            ):
                merged_line = LineString(
                    list(lines.pop(i).coords[::-1])[:-1] + list(merged_line.coords)
                )
                i = 0  # Restart the scan after merging
            else:
                i += 1
        merged_lines.append(merged_line)
        no_coordinates = [len(merged_lines[i].coords) for i in range(len(merged_lines))]
        max_index = np.argmax(no_coordinates)
        merged_lines = [merged_lines[max_index]]

    return gpd.GeoDataFrame(geometry=merged_lines, crs=gdf.crs)


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_osm_network")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # load default crs
    geo_crs = "EPSG:4326"
    distance_crs = "EPSG:3035"

    countries = snakemake.config["countries"]

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:
        build_network(
            snakemake.input,
            snakemake.output,
            geo_crs,
            distance_crs,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")
