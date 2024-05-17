# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import logging
import os

import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers import (
    configure_logging,
    set_scenario_config,
)
from shapely.geometry import LineString, Point
from shapely.ops import linemerge, split
from tqdm import tqdm
from _benchmark import memory_logger 

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
    if "na_rep" in kwargs:
        del kwargs["na_rep"]
    # if len(df) > 0:
    if not df.empty or not df.columns.empty:
        return df.to_csv(path, **kwargs, na_rep=NA_VALUES[0])
    else:
        with open(path, "w") as fp:
            pass


def line_endings_to_bus_conversion(lines):
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
        is_dclink_boundary_point = any(g_value["is_dclink_boundary_point"])

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
                    is_dclink_boundary_point, # check if new bus was formed of at least one DC link boundary point
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
        "is_dclink_boundary_point",
        "geometry",
    ]

    gdf_buses_clean = gpd.GeoDataFrame(buses_clean, columns=buses_clean_columns).set_crs(
        crs=buses.crs, inplace=True
    )

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


def get_converters(buses, lines):
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
                if g_value["is_dclink_boundary_point"].loc[id_0]:
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


# Note tolerance = 0.01 means around 700m
# TODO: the current tolerance is high to avoid an issue in the Nigeria case where line 565939360-1
#       seems to be interconnected to both ends, but at the eastern one, the node is actually not connected
#       another line seems to be exactly touching the node, but from the data point of view it only fly over it.
#       There may be the need to split a line in several segments in the case the line is within tolerance with
#       respect to a node


def merge_stations_lines_by_station_id_and_voltage(
    lines, buses, geo_crs, distance_crs, tol=5000
):
    """
    Function to merge close stations and adapt the line datasets to adhere to
    the merged dataset.
    """

    logger.info(
        "Stage 4a/5: Set substation ids with tolerance of %.2f km" % (tol / 1000)
    )

    # TODO pypsa-eur: Add this fix to pypsa-earth: Buses should not be clustered geographically if they are different 
    # bus types (AC != DC)
    buses_ac = buses[buses["dc"] == False].reset_index()
    buses_dc = buses[buses["dc"] == True].reset_index() 

    # set substation ids
    # set_substations_ids(buses, distance_crs, tol=tol)
    set_substations_ids(buses_ac, distance_crs, tol=tol)
    set_substations_ids(buses_dc, distance_crs, tol=tol)

    # Find boundary points of DC links
    # lines_dc_shape = lines[lines["dc"] == True].unary_union
    # lines_dc_bounds = lines_dc_shape.boundary
    # lines_dc_points = [p for p in lines_dc_bounds.geoms]
    lines_dc = lines[lines['dc'] == True].reset_index()
    lines_dc["adj_idx"] = range(0, len(lines_dc))

    # Initialize an empty adjacency matrix
    dc_adj_matrix = np.zeros((len(lines_dc), len(lines_dc)), dtype=int)

    # Fill the adjacency matrix
    for i in range(len(lines_dc)):
        for j in range(len(lines_dc)):
            if are_lines_connected(lines_dc.iloc[i], lines_dc.iloc[j]):
                dc_adj_matrix[i, j] = 1

    dc_paths = find_paths(dc_adj_matrix)

    all_dc_boundary_points = pd.Series()

    for path in dc_paths:
        bus_0_coors = lines_dc.iloc[path]["bus_0_coors"]
        bus_1_coors = lines_dc.iloc[path]["bus_1_coors"]

        # Create DataFrame containing all points within a path
        dc_points = pd.concat([bus_0_coors, bus_1_coors], ignore_index = True)

        # Determine the value counts of individual points. If it occurs more than 
        # once, it cannot be an end-point of a path
        bool_duplicates = dc_points.apply(lambda p: sum([are_almost_equal(p, s) for s in dc_points])) > 1
        
        # Drop all duplicates
        dc_boundary_points = dc_points[~bool_duplicates]

        if dc_boundary_points.empty:
            all_dc_boundary_points = dc_boundary_points
        else:
            all_dc_boundary_points = pd.concat([all_dc_boundary_points, dc_boundary_points], ignore_index = True)


    # TODO pypsa-eur: Add to pypsa-earth for all related entries on is_dclink_boundary_point
    # check for each entry in buses_dc whether it is included in lines_dc_points
    buses_ac["is_dclink_boundary_point"] = False
    buses_dc["is_dclink_boundary_point"] = buses_dc.geometry.apply(
        lambda p: any([p.within(l) for l in all_dc_boundary_points])
    )

    logger.info("Stage 4b/5: Merge substations with the same id")

    # merge buses with same station id and voltage
    if not buses.empty:
        buses_ac = merge_stations_same_station_id(buses_ac)
        buses_dc = merge_stations_same_station_id(buses_dc)
        buses_dc["bus_id"] = buses_ac["bus_id"].max() + buses_dc["bus_id"] + 1
        buses = pd.concat([buses_ac, buses_dc], ignore_index=True)
        set_substations_ids(buses, distance_crs, tol=tol)

    logger.info("Stage 4c/5: Specify the bus ids of the line endings")

    # set the bus ids to the line dataset
    lines, buses = set_lines_ids(lines, buses, distance_crs)

    # drop lines starting and ending in the same node
    lines.drop(lines[lines["bus0"] == lines["bus1"]].index, inplace=True)
    # update line endings
    lines = line_endings_to_bus_conversion(lines)

    # set substation_lv
    set_lv_substations(buses)

    logger.info("Stage 3d/4: Add converters to lines")

    # append fake converters
    # lines = pd.concat([lines, converters], ignore_index=True)

    # reset index
    lines.reset_index(drop=True, inplace=True)
    # if len(links) > 0:
    #     links.reset_index(drop=True, inplace=True)

    return lines, buses


def create_station_at_equal_bus_locations(
    lines, buses, geo_crs, distance_crs, tol=5000
):
    # V1. Create station_id at same bus location
    # - We saw that buses are not connected exactly at one point, they are
    #   usually connected to a substation "area" (analysed on maps)
    # - Create station_id at exactly the same location might therefore be not
    #   always correct
    # - Though as you can see below, it might be still sometime the case.
    #   Examples are **station 4** (2 lines with the same voltage connect at the
    #   same point) and **station 23** (4 lines with two different voltages connect
    #   at the same point)
    # TODO: Filter out the generator lines - defined as going from generator to
    #       the next station which is connected to a load. Excluding generator
    #       lines make probably sense because they are not transmission expansion
    #       relevant. For now we simplify and include generator lines.

    # If same location/geometry make station
    bus_all = buses

    # set substation ids
    set_substations_ids(buses, distance_crs, tol=tol)

    # set the bus ids to the line dataset
    lines, buses = set_lines_ids(lines, buses, distance_crs)

    # update line endings
    lines = line_endings_to_bus_conversion(lines)

    # For each station number with multiple buses make lowest voltage `substation_lv = TRUE`
    set_lv_substations(bus_all)

    # TRY: Keep only buses that are not duplicated & lv_substation = True
    # TODO: Check if this is necessary. What effect do duplicates have?
    bus_all = bus_all[bus_all["substation_lv"] == True]

    lines = connect_stations_same_station_id(lines, buses)

    return lines, buses


def _split_linestring_by_point(linestring, points):
    """
    Function to split a linestring geometry by multiple inner points.

    Parameters
    ----------
    lstring : LineString
        Linestring of the line to be split
    points : list
        List of points to split the linestring

    Return
    ------
    list_lines : list
        List of linestring to split the line
    """

    list_linestrings = [linestring]

    for p in points:
        # execute split to all lines and store results
        temp_list = [split(l, p) for l in list_linestrings]
        # nest all geometries
        list_linestrings = [lstring for tval in temp_list for lstring in tval.geoms]

    return list_linestrings


def fix_overpassing_lines(lines, buses, distance_crs, tol=1):
    """
    Function to avoid buses overpassing lines with no connection when the bus
    is within a given tolerance from the line.

    Parameters
    ----------
    lines : GeoDataFrame
        Geodataframe of lines
    buses : GeoDataFrame
        Geodataframe of substations
    tol : float
        Tolerance in meters of the distance between the substation and the line
        below which the line will be split
    """

    lines_to_add = []  # list of lines to be added
    lines_to_split = []  # list of lines that have been split

    lines_epsgmod = lines.to_crs(distance_crs)
    buses_epsgmod = buses.to_crs(distance_crs)

    # set tqdm options for substation ids
    tqdm_kwargs_substation_ids = dict(
        ascii=False,
        unit=" lines",
        total=lines.shape[0],
        desc="Verify lines overpassing nodes ",
    )

    for l in tqdm(lines.index, **tqdm_kwargs_substation_ids):
        # bus indices being within tolerance from the line
        bus_in_tol_epsg = buses_epsgmod[
            buses_epsgmod.geometry.distance(lines_epsgmod.geometry.loc[l]) <= tol
        ]

        # exclude endings of the lines
        bus_in_tol_epsg = bus_in_tol_epsg[
            (
                (
                    bus_in_tol_epsg.geometry.distance(
                        lines_epsgmod.geometry.loc[l].boundary.geoms[0]
                    )
                    > tol
                )
                | (
                    bus_in_tol_epsg.geometry.distance(
                        lines_epsgmod.geometry.loc[l].boundary.geoms[1]
                    )
                    > tol
                )
            )
        ]

        if not bus_in_tol_epsg.empty:
            # add index of line to split
            lines_to_split.append(l)

            buses_locs = buses.geometry.loc[bus_in_tol_epsg.index]

            # get new line geometries
            new_geometries = _split_linestring_by_point(lines.geometry[l], buses_locs)
            n_geoms = len(new_geometries)

            # create temporary copies of the line
            df_append = gpd.GeoDataFrame([lines.loc[l]] * n_geoms)
            # update geometries
            df_append["geometry"] = new_geometries
            # update name of the line
            df_append["line_id"] = [
                str(df_append["line_id"].iloc[0]) + f"_{id}" for id in range(n_geoms)
            ]

            lines_to_add.append(df_append)

    if not lines_to_add:
        return lines, buses

    df_to_add = gpd.GeoDataFrame(pd.concat(lines_to_add, ignore_index=True))
    df_to_add.set_crs(lines.crs, inplace=True)
    df_to_add.set_index(lines.index[-1] + df_to_add.index, inplace=True)

    # update length
    df_to_add["length"] = df_to_add.to_crs(distance_crs).geometry.length

    # update line endings
    df_to_add = line_endings_to_bus_conversion(df_to_add)

    # remove original lines
    lines.drop(lines_to_split, inplace=True)

    lines = gpd.GeoDataFrame(
        pd.concat([lines, df_to_add], ignore_index=True).reset_index(drop=True),
        crs=lines.crs,
    )

    return lines, buses


def add_buses_to_empty_countries(country_list, fp_country_shapes, buses):
    """
    Function to add a bus for countries missing substation data.
    """
    country_shapes = gpd.read_file(fp_country_shapes).set_index("name")["geometry"]
    bus_country_list = buses["country"].unique().tolist()

    # it may happen that bus_country_list contains entries not relevant as a country name (e.g. "not found")
    # difference can't give negative values; the following will return only relevant country names
    no_data_countries = list(set(country_list).difference(set(bus_country_list)))

    if len(no_data_countries) > 0:
        logger.info(
            f"No buses for the following countries: {no_data_countries}. Adding a node for everyone of them."
        )
        no_data_countries_shape = (
            country_shapes[country_shapes.index.isin(no_data_countries) == True]
            .reset_index()
            .to_crs(geo_crs)
        )
        length = len(no_data_countries)
        df = gpd.GeoDataFrame(
            {
                "voltage": [220000] * length,
                "country": no_data_countries_shape["name"],
                "x": no_data_countries_shape["geometry"].centroid.x,
                "y": no_data_countries_shape["geometry"].centroid.y,
                "bus_id": np.arange(len(buses) + 1, len(buses) + (length + 1), 1),
                "station_id": [np.nan] * length,
                # All lines for the countries with NA bus data are assumed to be AC
                "dc": [False] * length,
                "under_construction": [False] * length,
                "tag_area": [0.0] * length,
                "symbol": ["substation"] * length,
                "tag_substation": ["transmission"] * length,
                "geometry": no_data_countries_shape["geometry"].centroid,
                "substation_lv": [True] * length,
            },
            crs=geo_crs,
        ).astype(
            buses.dtypes.to_dict()
        )  # keep the same dtypes as buses
        buses = gpd.GeoDataFrame(
            pd.concat([buses, df], ignore_index=True).reset_index(drop=True),
            crs=buses.crs,
        )

        # update country list by buses dataframe
        bus_country_list = buses["country"].unique().tolist()

    non_allocated_countries = list(
        set(country_list).symmetric_difference(set(bus_country_list))
    )

    if len(non_allocated_countries) > 0:
        logger.error(
            f"There following countries could not be allocated properly: {non_allocated_countries}"
        )

    return buses


def build_network(
    inputs,
    outputs,
    build_osm_network_config,
    countries_config,
    geo_crs,
    distance_crs,
):  
    osm_clean_columns = {
        'substation': {
            'bus_id': 'object',
            'station_id': 'float',
            'voltage': 'float',
            'dc': 'bool',
            'symbol': 'object',
            'under_construction': 'bool',
            'tag_substation': 'str',
            'tag_area': 'str',
            'lon': 'float',
            'lat': 'float',
            'country': 'str',
            'geometry': 'object',
            'tag_source': 'str',
        },
        'line': {
            'line_id': 'object',
            'bus0': 'object',
            'bus1': 'object',
            'voltage': 'float',
            'circuits': 'float',
            'length': 'float',
            'underground': 'bool',
            'under_construction': 'bool',
            'tag_type': 'str',
            'tag_frequency': 'float',
            'dc': 'bool',
            'country': 'object',
            'geometry': 'object',
        }
    }

    logger.info("Stage 1/5: Read input data")
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

    lines = line_endings_to_bus_conversion(lines)

    logger.info("Stage 2/5: AC and DC network: enabled")

    # Address the overpassing line issue Step 3/5
    if build_osm_network_config.get("split_overpassing_lines", False):
        tol = build_osm_network_config.get("overpassing_lines_tolerance", 1)
        logger.info("Stage 3/5: Avoid nodes overpassing lines: enabled with tolerance")

        lines, buses = fix_overpassing_lines(lines, buses, distance_crs, tol=tol)
    else:
        logger.info("Stage 3/5: Avoid nodes overpassing lines: disabled")
    
    # Add bus to countries with no buses
    buses = add_buses_to_empty_countries(countries_config, inputs.country_shapes, buses)

    # METHOD to merge buses with same voltage and within tolerance Step 4/5
    tol = build_osm_network_config.get("group_tolerance_buses", 5000)
    logger.info(
        f"Stage 4/5: Aggregate close substations: enabled with tolerance {tol} m"
    )
    lines, buses = merge_stations_lines_by_station_id_and_voltage(
        lines, buses, geo_crs, distance_crs, tol=tol
    )

    logger.info("Stage 5/5: Add augmented substation to country with no data")

    # Recalculate lengths of lines
    utm = lines.estimate_utm_crs(datum_name = "WGS 84")
    lines["length"] = lines.to_crs(utm).length

    # get transformers: modelled as lines connecting buses with different voltage
    transformers = get_transformers(buses, lines)

    # get converters: currently modelled as links connecting buses with different polarity
    converters = get_converters(buses, lines)

    logger.info("Save outputs")

    # create clean directory if not already exist
    if not os.path.exists(outputs["lines"]):
        os.makedirs(os.path.dirname(outputs["lines"]), exist_ok=True)


    ### Convert output to pypsa-eur friendly format
    # Rename "substation" in buses["symbol"] to "Substation"
    buses["symbol"] = buses["symbol"].replace({"substation": "Substation"})

    # Drop unncessary index column and set respective element ids as index
    lines.set_index("line_id", inplace=True)
    converters.set_index("converter_id", inplace=True)
    transformers.set_index("transformer_id", inplace=True)
    buses.set_index("bus_id", inplace=True)


    # Convert voltages from V to kV
    lines["voltage"] = lines["voltage"] / 1000
    transformers["voltage_bus0"], transformers["voltage_bus1"] = transformers["voltage_bus0"] / 1000, \
        transformers["voltage_bus1"] / 1000
    buses["voltage"] = buses["voltage"] / 1000

    # Convert 'true' and 'false' to 't' and 'f'    
    lines = lines.replace({True: "t", False: "f"})
    converters = converters.replace({True: "t", False: "f"})
    buses = buses.replace({True: "t", False: "f"})
    
    # Change column orders
    cols_lines = ["bus0", "bus1", "voltage", "circuits", "length", "underground", "under_construction", "geometry",
                  "tag_type", "tag_frequency", "country", "bounds", 
                  "bus_0_coors", "bus_1_coors", "bus0_lon", "bus0_lat", "bus1_lon", "bus1_lat"]
    
    cols_lines_csv = ["bus0", "bus1", "voltage", "circuits", "tag_frequency", "length", "underground", "under_construction", "geometry"]
    lines_csv = lines[cols_lines_csv]
    lines = lines[cols_lines]
    


    to_csv_nafix(lines_csv, outputs["lines"], quotechar="'")  # Generate CSV
    to_csv_nafix(converters, outputs["converters"], quotechar="'")  # Generate CSV
    to_csv_nafix(transformers, outputs["transformers"], quotechar="'")  # Generate CSV

    colstodrop = ["bounds", "bus_0_coors", "bus_1_coors"]

    # Export to GeoJSON for quick validations
    save_to_geojson(gpd.GeoDataFrame(lines.drop(columns = colstodrop), geometry = "geometry", crs = geo_crs), outputs["lines_geojson"])
    save_to_geojson(gpd.GeoDataFrame(converters, geometry = "geometry", crs = geo_crs), outputs["converters_geojson"])
    save_to_geojson(gpd.GeoDataFrame(transformers.drop(columns = colstodrop), geometry = "geometry", crs = geo_crs), outputs["transformers_geojson"])

    # create clean directory if not already exist
    if not os.path.exists(outputs["substations"]):
        os.makedirs(os.path.dirname(outputs["substations"]), exist_ok=True)
    # Generate CSV
    to_csv_nafix(buses, outputs["substations"], quotechar="'")
    save_to_geojson(gpd.GeoDataFrame(buses, geometry = "geometry", crs = geo_crs), outputs["substations_geojson"])

    return None


# Function to check if two lines are connected
def are_lines_connected(line1, line2):
    # return (line1['geometry'].touches(line2['geometry']))
    return (
        are_almost_equal(line1["bus_0_coors"], line2["bus_0_coors"]),
        are_almost_equal(line1["bus_0_coors"], line2["bus_1_coors"]),
        are_almost_equal(line1["bus_1_coors"], line2["bus_0_coors"]),
        are_almost_equal(line1["bus_1_coors"], line2["bus_1_coors"])
        )


def _dfs(adj_matrix, visited, current_vertex, path):
    visited[current_vertex] = True
    path.append(current_vertex)
    for neighbor in range(len(adj_matrix)):
        if adj_matrix[current_vertex][neighbor] == 1 and not visited[neighbor]:
            _dfs(adj_matrix, visited, neighbor, path)
    return path


# Returns all connected paths as a vector
def find_paths(adj_matrix):
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

    build_osm_network = snakemake.config["build_osm_network"]
    countries = snakemake.config["countries"]

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:
        build_network(
        snakemake.input,
        snakemake.output,
        build_osm_network,
        countries,
        geo_crs,
        distance_crs,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")