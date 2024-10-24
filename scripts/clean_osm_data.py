# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script is used to clean OpenStreetMap (OSM) data for creating a PyPSA-Eur
ready network.

The script performs various cleaning operations on the OSM data, including:
- Cleaning voltage, circuits, cables, wires, and frequency columns
- Splitting semicolon-separated cells into new rows
- Distributing values to circuits based on the number of splits
- Adding line endings to substations based on line data
"""

import itertools
import json
import logging
import os
import re

import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from shapely.algorithms.polylabel import polylabel
from shapely.geometry import LineString, MultiLineString, Point, Polygon
from shapely.ops import linemerge, unary_union

logger = logging.getLogger(__name__)


GEO_CRS = "EPSG:4326"
DISTANCE_CRS = "EPSG:3035"
BUS_TOL = (
    500  # unit: meters, default 5000 - Buses within this distance are grouped together
)


def _create_linestring(row):
    """
    Create a LineString object from the given row.

    Args:
        row (dict): A dictionary containing the row data.

    Returns:
        LineString: A LineString object representing the geometry.
    """
    coords = [(coord["lon"], coord["lat"]) for coord in row["geometry"]]
    return LineString(coords)


def _create_polygon(row):
    """
    Create a Shapely Polygon from a list of coordinate dictionaries.

    Parameters:
        coords (list): List of dictionaries with 'lat' and 'lon' keys
        representing coordinates.

    Returns:
        shapely.geometry.Polygon: The constructed polygon object.
    """
    # Extract coordinates as tuples
    point_coords = [(coord["lon"], coord["lat"]) for coord in row["geometry"]]

    # Ensure closure by repeating the first coordinate as the last coordinate
    if point_coords[0] != point_coords[-1]:
        point_coords.append(point_coords[0])

    # Create Polygon object
    polygon = Polygon(point_coords)

    return polygon


def _clean_voltage(column):
    """
    Function to clean the raw voltage column: manual fixing and drop nan values

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
    logger.info("Cleaning voltages.")
    column = column.copy()

    column = (
        column.astype(str)
        .str.lower()
        .str.replace("400/220/110 kV'", "400000;220000;110000")
        .str.replace("400/220/110/20_kv", "400000;220000;110000;20000")
        .str.replace("2x25000", "25000;25000")
        .str.replace("é", ";")
    )

    column = (
        column.astype(str)
        .str.lower()
        .str.replace("(temp 150000)", "")
        .str.replace("low", "1000")
        .str.replace("minor", "1000")
        .str.replace("medium", "33000")
        .str.replace("med", "33000")
        .str.replace("m", "33000")
        .str.replace("high", "150000")
        .str.replace("23000-109000", "109000")
        .str.replace("380000>220000", "380000;220000")
        .str.replace(":", ";")
        .str.replace("<", ";")
        .str.replace(",", ";")
        .str.replace("kv", "000")
        .str.replace("kva", "000")
        .str.replace("/", ";")
        .str.replace("nan", "")
        .str.replace("<na>", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r"[^0-9;]", "", str(x)))

    column.dropna(inplace=True)
    return column


def _clean_circuits(column):
    """
    Function to clean the raw circuits column: manual fixing and drop nan
    values

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
    logger.info("Cleaning circuits.")
    column = column.copy()
    column = (
        column.astype(str)
        .str.replace("partial", "")
        .str.replace("1operator=RTE operator:wikidata=Q2178795", "")
        .str.lower()
        .str.replace("1,5", "3")
        .str.replace("1/3", "1")
        .str.replace("<na>", "")
        .str.replace("nan", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r"[^0-9;]", "", x))

    column.dropna(inplace=True)
    return column.astype(str)


def _clean_cables(column):
    """
    Function to clean the raw cables column: manual fixing and drop nan values

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
    logger.info("Cleaning cables.")
    column = column.copy()
    column = (
        column.astype(str)
        .str.lower()
        .str.replace("1/3", "1")
        .str.replace("3x2;2", "3")
        .str.replace("<na>", "")
        .str.replace("nan", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r"[^0-9;]", "", x))

    column.dropna(inplace=True)
    return column.astype(str)


def _clean_wires(column):
    """
    Function to clean the raw wires column: manual fixing and drop nan values

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
    logger.info("Cleaning wires.")
    column = column.copy()
    column = (
        column.astype(str)
        .str.lower()
        .str.replace("?", "")
        .str.replace("trzyprzewodowe", "3")
        .str.replace("pojedyńcze", "1")
        .str.replace("single", "1")
        .str.replace("double", "2")
        .str.replace("triple", "3")
        .str.replace("quad", "4")
        .str.replace("fivefold", "5")
        .str.replace("yes", "3")
        .str.replace("1/3", "1")
        .str.replace("3x2;2", "3")
        .str.replace("_", "")
        .str.replace("<na>", "")
        .str.replace("nan", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r"[^0-9;]", "", x))

    column.dropna(inplace=True)
    return column.astype(str)


def _check_voltage(voltage, list_voltages):
    """
    Check if the given voltage is present in the list of allowed voltages.

    Parameters:
    voltage (str): The voltage to check.
    list_voltages (list): A list of allowed voltages.

    Returns:
    bool: True if the voltage is present in the list of allowed voltages,
    False otherwise.
    """
    voltages = voltage.split(";")
    for v in voltages:
        if v in list_voltages:
            return True
    return False


def _clean_frequency(column):
    """
    Function to clean the raw frequency column: manual fixing and drop nan
    values

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
    logger.info("Cleaning frequencies.")
    column = column.copy()
    column = (
        column.astype(str)
        .str.lower()
        .str.replace("16.67", "16.7")
        .str.replace("16,7", "16.7")
        .str.replace("?", "")
        .str.replace("hz", "")
        .str.replace(" ", "")
        .str.replace("<NA>", "")
        .str.replace("nan", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r"[^0-9;.]", "", x))

    column.dropna(inplace=True)
    return column.astype(str)


def _clean_rating(column):
    """
    Function to clean and sum the rating columns:

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
    logger.info("Cleaning ratings.")
    column = column.copy()
    column = column.astype(str).str.replace("MW", "")

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r"[^0-9;]", "", x))

    # Sum up all ratings if there are multiple entries
    column = column.str.split(";").apply(lambda x: sum([int(i) for i in x]))

    column.dropna(inplace=True)
    return column.astype(str)


def _split_cells(df, cols=["voltage"]):
    """
    Split semicolon separated cells i.e. [66000;220000] and create new
    identical rows.

    Parameters
    ----------
    df : dataframe
        Dataframe under analysis
    cols : list
        List of target columns over which to perform the analysis

    Example
    -------
    Original data:
    row 1: '66000;220000', '50'

    After applying split_cells():
    row 1, '66000', '50', 2
    row 2, '220000', '50', 2
    """
    if df.empty:
        return df

    # Create a dictionary to store the suffix count for each original ID
    suffix_counts = {}
    # Create a dictionary to store the number of splits associated with each
    # original ID
    num_splits = {}

    # Split cells and create new rows
    x = df.assign(**{col: df[col].str.split(";") for col in cols})
    x = x.explode(cols, ignore_index=True)

    # Count the number of splits associated with each original ID
    num_splits = x.groupby("id").size().to_dict()

    # Update the 'split_elements' column
    x["split_elements"] = x["id"].map(num_splits)

    # Function to generate the new ID with suffix and update the number of
    # splits
    def generate_new_id(row):
        original_id = row["id"]
        if row["split_elements"] == 1:
            return original_id
        else:
            suffix_counts[original_id] = suffix_counts.get(original_id, 0) + 1
            return f"{original_id}-{suffix_counts[original_id]}"

    # Update the ID column with the new IDs
    x["id"] = x.apply(generate_new_id, axis=1)

    return x


def _distribute_to_circuits(row):
    """
    Distributes the number of circuits or cables to individual circuits based
    on the given row data.

    Parameters:
    - row: A dictionary representing a row of data containing information about
      circuits and cables.

    Returns:
    - single_circuit: The number of circuits to be assigned to each individual
      circuit.
    """
    if row["circuits"] != "":
        circuits = int(row["circuits"])
    else:
        cables = int(row["cables"])
        circuits = cables / 3

    single_circuit = int(max(1, np.floor_divide(circuits, row["split_elements"])))
    single_circuit = str(single_circuit)

    return single_circuit


def _import_lines_and_cables(path_lines):
    """
    Import lines and cables from the given input paths.

    Parameters:
    - path_lines (dict): A dictionary containing the input paths for lines and
      cables data.

    Returns:
    - df_lines (DataFrame): A DataFrame containing the imported lines and
      cables data.
    """
    columns = [
        "id",
        "bounds",
        "nodes",
        "geometry",
        "country",
        "power",
        "cables",
        "circuits",
        "frequency",
        "voltage",
        "wires",
    ]
    df_lines = pd.DataFrame(columns=columns)

    logger.info("Importing lines and cables")
    for key in path_lines:
        logger.info(f"Processing {key}...")
        for idx, ip in enumerate(path_lines[key]):
            if (
                os.path.exists(ip) and os.path.getsize(ip) > 400
            ):  # unpopulated OSM json is about 51 bytes
                country = os.path.basename(os.path.dirname(path_lines[key][idx]))

                logger.info(
                    f" - Importing {key} {str(idx+1).zfill(2)}/{str(len(path_lines[key])).zfill(2)}: {ip}"
                )
                with open(ip, "r") as f:
                    data = json.load(f)

                df = pd.DataFrame(data["elements"])
                df["id"] = df["id"].astype(str)
                df["country"] = country

                col_tags = [
                    "power",
                    "cables",
                    "circuits",
                    "frequency",
                    "voltage",
                    "wires",
                ]

                tags = pd.json_normalize(df["tags"]).map(
                    lambda x: str(x) if pd.notnull(x) else x
                )

                for ct in col_tags:
                    if ct not in tags.columns:
                        tags[ct] = pd.NA

                tags = tags.loc[:, col_tags]

                df = pd.concat([df, tags], axis="columns")
                df.drop(columns=["type", "tags"], inplace=True)

                df_lines = pd.concat([df_lines, df], axis="rows")

            else:
                logger.info(
                    f" - Skipping {key} {str(idx+1).zfill(2)}/{str(len(path_lines[key])).zfill(2)} (empty): {ip}"
                )
                continue
        logger.info("---")

    # Append prefix "way/"
    df_lines["id"] = "way/" + df_lines["id"]

    return df_lines


def _import_routes_relation(path_relation):
    """ """
    columns = [
        "id",
        "bounds",
        "nodes",
        "geometry",
        "country",
        "circuits",
        "cables",
        "frequency",
        "voltage",
    ]
    df_relation = pd.DataFrame(columns=columns)

    logger.info("Importing power route relations (lines, cables, links)")
    for key in path_relation:
        logger.info(f"Processing {key}...")
        for idx, ip in enumerate(path_relation[key]):
            if (
                os.path.exists(ip) and os.path.getsize(ip) > 400
            ):  # unpopulated OSM json is about 51 bytes
                country = os.path.basename(os.path.dirname(path_relation[key][idx]))

                logger.info(
                    f" - Importing {key} {str(idx+1).zfill(2)}/{str(len(path_relation[key])).zfill(2)}: {ip}"
                )
                with open(ip, "r") as f:
                    data = json.load(f)

                df = pd.DataFrame(data["elements"])
                df["id"] = df["id"].astype(str)
                df["id"] = df["id"].apply(lambda x: (f"relation/{x}"))
                df["country"] = country

                col_tags = [
                    "circuits",
                    "cables",
                    "frequency",
                    "voltage",
                    "rating",
                ]

                tags = pd.json_normalize(df["tags"]).map(
                    lambda x: str(x) if pd.notnull(x) else x
                )

                for ct in col_tags:
                    if ct not in tags.columns:
                        tags[ct] = pd.NA

                tags = tags.loc[:, col_tags]

                df = pd.concat([df, tags], axis="columns")
                df.drop(columns=["type", "tags"], inplace=True)

                df_relation = pd.concat([df_relation, df], axis="rows")

            else:
                logger.info(
                    f" - Skipping {key} {str(idx+1).zfill(2)}/{str(len(path_relation[key])).zfill(2)} (empty): {ip}"
                )
                continue

    return df_relation


def _create_single_link(row):
    """
    Create a single link from multiple rows within a OSM link relation.

    Parameters:
    - row: A row of OSM data containing information about the link.

    Returns:
    - single_link: A single LineString representing the link.

    This function takes a row of OSM data and extracts the relevant information
    to create a single link. It filters out elements (substations, electrodes)
    with invalid roles and finds the longest link based on its endpoints.
    If the longest link is a MultiLineString, it extracts the longest
    linestring from it. The resulting single link is returned.
    """
    valid_roles = ["line", "cable"]
    df = pd.json_normalize(row["members"])
    df = df[df["role"].isin(valid_roles)]
    df.loc[:, "geometry"] = df.apply(_create_linestring, axis=1)
    df.loc[:, "length"] = df["geometry"].apply(lambda x: x.length)

    list_endpoints = []
    for idx, row in df.iterrows():
        tuple = sorted([row["geometry"].coords[0], row["geometry"].coords[-1]])
        # round tuple to 3 decimals
        tuple = (
            round(tuple[0][0], 3),
            round(tuple[0][1], 3),
            round(tuple[1][0], 3),
            round(tuple[1][1], 3),
        )
        list_endpoints.append(tuple)

    df.loc[:, "endpoints"] = list_endpoints
    df_longest = df.loc[df.groupby("endpoints")["length"].idxmin()]

    single_link = linemerge(df_longest["geometry"].values.tolist())

    # If the longest component is a MultiLineString, extract the longest linestring from it
    if isinstance(single_link, MultiLineString):
        # Find connected components
        components = list(single_link.geoms)

        # Find the longest connected linestring
        single_link = max(components, key=lambda x: x.length)

    return single_link


def _extract_members(row):
    df = pd.json_normalize(row["members"])
    df["ref"] = df["ref"].astype(str)
    df["ways"] = "way/" + df["ref"]
    member_ids = df["ways"].values.tolist()

    return member_ids


def _create_line(row):
    """
    Create a line from multiple rows. Drops closed geometries (substations).

    Parameters:
    - row: A row of OSM data containing information about the relation/line.

    Returns:
    - line: LineString/MultiLineString representing the relation/line.
    """
    df = pd.json_normalize(row["members"])
    df["ref"] = df["ref"].astype(str)
    df["ways"] = "way/" + df["ref"]
    # Drop NAs
    df = df.dropna(subset=["geometry"])
    df.loc[:, "geometry"] = df.apply(_create_linestring, axis=1)
    # Drop closed geometries (substations)
    closed_geom = df["geometry"].apply(lambda x: x.is_closed)

    line = linemerge(df[~closed_geom]["geometry"].values.tolist())
    members = df[~closed_geom]["ways"].values.tolist()

    return line, members


def _drop_duplicate_lines(df_lines):
    """
    Drop duplicate lines from the given dataframe. Duplicates are usually lines
    cross-border lines or slightly outside the country border of focus.

    Parameters:
    - df_lines (pandas.DataFrame): The dataframe containing lines data.

    Returns:
    - df_lines (pandas.DataFrame): The dataframe with duplicate lines removed
      and cleaned data.

    This function drops duplicate lines from the given dataframe based on the
    'id' column. It groups the duplicate rows by 'id' and aggregates the
    'country' column to a string split by semicolon, as they appear in multiple
    country datasets. One example of the duplicates is kept, accordingly.
    Finally, the updated dataframe without multiple duplicates is returned.
    """
    logger.info("Dropping duplicate lines.")
    duplicate_rows = df_lines[df_lines.duplicated(subset=["id"], keep=False)].copy()

    # Group rows by id and aggregate the country column to a string split by semicolon
    grouped_duplicates = (
        duplicate_rows.groupby("id")["country"].agg(lambda x: ";".join(x)).reset_index()
    )
    duplicate_rows.drop_duplicates(subset="id", inplace=True)
    duplicate_rows.drop(columns=["country"], inplace=True)
    duplicate_rows = duplicate_rows.join(
        grouped_duplicates.set_index("id"), on="id", how="left"
    )

    len_before = len(df_lines)
    # Drop duplicates and update the df_lines dataframe with the cleaned data
    df_lines = df_lines[~df_lines["id"].isin(duplicate_rows["id"])]
    df_lines = pd.concat([df_lines, duplicate_rows], axis="rows")
    len_after = len(df_lines)

    logger.info(
        f"Dropped {len_before - len_after} duplicate elements. "
        + f"Keeping {len_after} elements."
    )

    return df_lines


def _filter_by_voltage(df, min_voltage=220000):
    """
    Filter rows in the DataFrame based on the voltage in V.

    Parameters:
    - df (pandas.DataFrame): The DataFrame containing the substations or lines data.
    - min_voltage (int, optional): The minimum voltage value to filter the
      rows. Defaults to 220000 [unit: V].

    Returns:
    - filtered df (pandas.DataFrame): The filtered DataFrame containing
      the lines or substations above min_voltage.
    - list_voltages (list): A list of unique voltage values above min_voltage.
      The type of the list elements is string.
    """
    if df.empty:
        return df, []

    logger.info(
        f"Filtering dataframe by voltage. Only keeping rows above and including {min_voltage} V."
    )
    list_voltages = df["voltage"].str.split(";").explode().unique().astype(str)
    # Keep numeric strings
    list_voltages = list_voltages[np.vectorize(str.isnumeric)(list_voltages)]
    list_voltages = list_voltages.astype(int)
    list_voltages = list_voltages[list_voltages >= int(min_voltage)]
    list_voltages = list_voltages.astype(str)

    bool_voltages = df["voltage"].apply(_check_voltage, list_voltages=list_voltages)
    len_before = len(df)
    df = df[bool_voltages]
    len_after = len(df)
    logger.info(
        f"Dropped {len_before - len_after} elements with voltage below {min_voltage}. "
        + f"Keeping {len_after} elements."
    )

    return df, list_voltages


def _clean_substations(df_substations, list_voltages):
    """
    Clean the substation data by performing the following steps:
    - Split cells in the dataframe.
    - Filter substation data based on specified voltages.
    - Update the frequency values based on the split count.
    - Split cells in the 'frequency' column.
    - Set remaining invalid frequency values that are not in ['0', '50']
      to '50'.

    Parameters:
    - df_substations (pandas.DataFrame): The input dataframe containing
      substation data.
    - list_voltages (list): A list of voltages above min_voltage to filter the
    substation data.

    Returns:
    - df_substations (pandas.DataFrame): The cleaned substation dataframe.
    """
    df_substations = df_substations.copy()

    df_substations = _split_cells(df_substations)

    bool_voltages = df_substations["voltage"].apply(
        _check_voltage, list_voltages=list_voltages
    )
    df_substations = df_substations[bool_voltages]
    df_substations.loc[:, "split_count"] = df_substations["id"].apply(
        lambda x: x.split("-")[1] if "-" in x else "0"
    )
    df_substations.loc[:, "split_count"] = df_substations["split_count"].astype(int)

    bool_split = df_substations["split_elements"] > 1
    bool_frequency_len = (
        df_substations["frequency"].apply(lambda x: len(x.split(";")))
        == df_substations["split_elements"]
    )

    op_freq = lambda row: row["frequency"].split(";")[row["split_count"] - 1]

    df_substations.loc[bool_frequency_len & bool_split, "frequency"] = (
        df_substations.loc[bool_frequency_len & bool_split,].apply(op_freq, axis=1)
    )

    df_substations = _split_cells(df_substations, cols=["frequency"])
    bool_invalid_frequency = df_substations["frequency"].apply(
        lambda x: x not in ["50", "0"]
    )
    df_substations.loc[bool_invalid_frequency, "frequency"] = "50"

    return df_substations


def _clean_lines(df_lines, list_voltages):
    """
    Cleans and processes the `df_lines` DataFrame heuristically based on the
    information available per respective line and cable. Further checks to
    ensure data consistency and completeness.

    Parameters
    ----------
    df_lines : pandas.DataFrame
        The input DataFrame containing line information with columns such as
        'voltage', 'circuits', 'frequency', 'cables', 'split_elements', 'id',
        etc.
    list_voltages : list
        A list of unique voltage values above a certain threshold. (type: str)

    Returns
    -------
    df_lines : pandas.DataFrame
        The cleaned DataFrame with updated columns 'circuits', 'frequency', and
        'cleaned' to reflect the applied transformations.

    Description
    -----------
    This function performs the following operations:

    - Initializes a 'cleaned' column with False, step-wise updates to True
       following the respective cleaning step.
    - Splits the voltage cells in the DataFrame at semicolons using a helper
       function `_split_cells`.
    - Filters the DataFrame to only include rows with valid voltages.
    - Sets circuits of remaining lines without any applicable heuristic equal
      to 1.

    The function ensures that the resulting DataFrame has consistent and
    complete information for further processing or analysis while maintaining
    the data of the original OSM data set wherever possible.
    """
    logger.info("Cleaning lines and determining circuits.")
    # Initiate boolean with False, only set to true if all cleaning steps are
    # passed
    df_lines = df_lines.copy()
    df_lines["cleaned"] = False

    df_lines["voltage_original"] = df_lines["voltage"]
    df_lines["circuits_original"] = df_lines["circuits"]

    df_lines = _split_cells(df_lines)
    bool_voltages = df_lines["voltage"].apply(
        _check_voltage, list_voltages=list_voltages
    )
    df_lines = df_lines[bool_voltages]

    bool_ac = df_lines["frequency"] != "0"
    bool_dc = ~bool_ac
    valid_frequency = ["50", "0"]
    bool_invalid_frequency = df_lines["frequency"].apply(
        lambda x: x not in valid_frequency
    )

    bool_noinfo = (df_lines["cables"] == "") & (df_lines["circuits"] == "")
    # Fill in all values where cables info and circuits does not exist. Assuming 1 circuit
    df_lines.loc[bool_noinfo, "circuits"] = "1"
    df_lines.loc[bool_noinfo & bool_invalid_frequency, "frequency"] = "50"
    df_lines.loc[bool_noinfo, "cleaned"] = True

    # Fill in all values where cables info exists and split_elements == 1
    bool_cables_ac = (
        (df_lines["cables"] != "")
        & (df_lines["split_elements"] == 1)
        & (df_lines["cables"] != "0")
        & (df_lines["cables"].apply(lambda x: len(x.split(";")) == 1))
        & (df_lines["circuits"] == "")
        & (df_lines["cleaned"] == False)
        & bool_ac
    )

    df_lines.loc[bool_cables_ac, "circuits"] = df_lines.loc[
        bool_cables_ac, "cables"
    ].apply(lambda x: str(int(max(1, np.floor_divide(int(x), 3)))))

    df_lines.loc[bool_cables_ac, "frequency"] = "50"
    df_lines.loc[bool_cables_ac, "cleaned"] = True

    bool_cables_dc = (
        (df_lines["cables"] != "")
        & (df_lines["split_elements"] == 1)
        & (df_lines["cables"] != "0")
        & (df_lines["cables"].apply(lambda x: len(x.split(";")) == 1))
        & (df_lines["circuits"] == "")
        & (df_lines["cleaned"] == False)
        & bool_dc
    )

    df_lines.loc[bool_cables_dc, "circuits"] = df_lines.loc[
        bool_cables_dc, "cables"
    ].apply(lambda x: str(int(max(1, np.floor_divide(int(x), 2)))))

    df_lines.loc[bool_cables_dc, "frequency"] = "0"
    df_lines.loc[bool_cables_dc, "cleaned"] = True

    # Fill in all values where circuits info exists and split_elements == 1
    bool_lines = (
        (df_lines["circuits"] != "")
        & (df_lines["split_elements"] == 1)
        & (df_lines["circuits"] != "0")
        & (df_lines["circuits"].apply(lambda x: len(x.split(";")) == 1))
        & (df_lines["cleaned"] == False)
    )

    df_lines.loc[bool_lines & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_lines & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_lines, "cleaned"] = True

    # Clean those values where number of voltages split by semicolon is larger
    # than no cables or no circuits
    bool_cables = (
        (df_lines["voltage_original"].apply(lambda x: len(x.split(";")) > 1))
        & (df_lines["cables"].apply(lambda x: len(x.split(";")) == 1))
        & (df_lines["circuits"].apply(lambda x: len(x.split(";")) == 1))
        & (df_lines["cleaned"] == False)
    )

    df_lines.loc[bool_cables, "circuits"] = df_lines[bool_cables].apply(
        _distribute_to_circuits, axis=1
    )
    df_lines.loc[bool_cables & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_cables & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_cables, "cleaned"] = True

    # Clean those values where multiple circuit values are present, divided by
    # semicolon
    has_multiple_circuits = df_lines["circuits"].apply(lambda x: len(x.split(";")) > 1)
    circuits_match_split_elements = df_lines.apply(
        lambda row: len(row["circuits"].split(";")) == row["split_elements"],
        axis=1,
    )
    is_not_cleaned = df_lines["cleaned"] == False
    bool_cables = has_multiple_circuits & circuits_match_split_elements & is_not_cleaned

    df_lines.loc[bool_cables, "circuits"] = df_lines.loc[bool_cables].apply(
        lambda row: str(row["circuits"].split(";")[int(row["id"].split("-")[-1]) - 1]),
        axis=1,
    )

    df_lines.loc[bool_cables & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_cables & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_cables, "cleaned"] = True

    # Clean those values where multiple cables values are present, divided by
    # semicolon
    has_multiple_cables = df_lines["cables"].apply(lambda x: len(x.split(";")) > 1)
    cables_match_split_elements = df_lines.apply(
        lambda row: len(row["cables"].split(";")) == row["split_elements"],
        axis=1,
    )
    is_not_cleaned = df_lines["cleaned"] == False
    bool_cables = has_multiple_cables & cables_match_split_elements & is_not_cleaned

    df_lines.loc[bool_cables, "circuits"] = df_lines.loc[bool_cables].apply(
        lambda row: str(
            max(
                1,
                np.floor_divide(
                    int(row["cables"].split(";")[int(row["id"].split("-")[-1]) - 1]), 3
                ),
            )
        ),
        axis=1,
    )

    df_lines.loc[bool_cables & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_cables & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_cables, "cleaned"] = True

    # All remaining lines to circuits == 1
    bool_leftover = df_lines["cleaned"] == False
    if sum(bool_leftover) > 0:
        str_id = "; ".join(str(id) for id in df_lines.loc[bool_leftover, "id"])
        logger.info(f"Setting circuits of remaining {sum(bool_leftover)} lines to 1...")
        logger.info(f"Lines affected: {str_id}")

    df_lines.loc[bool_leftover, "circuits"] = "1"
    df_lines.loc[bool_leftover & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_leftover & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_leftover, "cleaned"] = True

    return df_lines


def _create_substations_geometry(df_substations):
    """
    Creates geometries.

    Parameters:
    df_substations (DataFrame): The input DataFrame containing the substations
    data.

    Returns:
    df_substations (DataFrame): A new DataFrame with the
    polygons ["polygon"] of the substations geometries.
    """
    logger.info("Creating substations geometry.")
    df_substations = df_substations.copy()

    # Create PoI from geometries and keep the original polygons
    df_substations.loc[:, "polygon"] = df_substations["geometry"]

    return df_substations


def _create_substations_poi(df_substations, tol=BUS_TOL / 2):
    """
    Creates Pole of Inaccessibility (PoI) from geometries and keeps the original polygons.

    Parameters:
    df_substations (DataFrame): The input DataFrame containing the substations
    data.

    Returns:
    df_substations (DataFrame): A new DataFrame with the PoI ["geometry"]
    and polygons ["polygon"] of the substations geometries.
    """
    logger.info("Creating substations geometry.")
    df_substations = df_substations.copy()

    df_substations.loc[:, "geometry"] = df_substations["polygon"].apply(
        lambda polygon: polylabel(polygon, tol)
    )

    df_substations.loc[:, "lon"] = df_substations["geometry"].apply(lambda x: x.x)
    df_substations.loc[:, "lat"] = df_substations["geometry"].apply(lambda x: x.y)

    return df_substations


def _create_lines_geometry(df_lines):
    """
    Create line geometry for the given DataFrame of lines.

    Parameters:
    - df_lines (pandas.DataFrame): DataFrame containing lines data.

    Returns:
    - df_lines (pandas.DataFrame): DataFrame with transformed 'geometry'
      column (type: shapely LineString).

    Notes:
    - This function transforms 'geometry' column in the input DataFrame by
      applying the '_create_linestring' function to each row.
    - It then drops rows where the geometry has equal start and end points,
      as these are usually not lines but outlines of areas.
    """
    logger.info("Creating lines geometry.")
    df_lines = df_lines.copy()
    df_lines.loc[:, "geometry"] = df_lines.apply(_create_linestring, axis=1)

    bool_circle = df_lines["geometry"].apply(lambda x: x.coords[0] == x.coords[-1])
    df_lines = df_lines[~bool_circle]

    return df_lines


def _add_bus_poi_to_line(linestring, point):
    """
    Adds the PoI of a substation to a linestring by extending the
    linestring with a new segment.

    Parameters:
    linestring (LineString): The original linestring to extend.
    point (Point): The PoI of the bus.

    Returns:
    merged (LineString): The extended linestring with the new segment.
    """
    start = linestring.coords[0]
    end = linestring.coords[-1]

    dist_to_start = point.distance(Point(start))
    dist_to_end = point.distance(Point(end))

    if dist_to_start < dist_to_end:
        new_segment = LineString([point.coords[0], start])
    else:
        new_segment = LineString([point.coords[0], end])

    merged = linemerge([linestring, new_segment])

    return merged


def _finalise_substations(df_substations):
    """
    Finalises the substations column types.

    Args:
        df_substations (pandas.DataFrame): The input DataFrame
        containing substations data.

    Returns:
        df_substations (pandas.DataFrame(): The DataFrame with finalised column
        types and transformed data.
    """
    logger.info("Finalising substations column types.")
    df_substations = df_substations.copy()
    # rename columns
    df_substations.rename(
        columns={
            "id": "bus_id",
            "power": "symbol",
            "substation": "tag_substation",
        },
        inplace=True,
    )

    # Initiate new columns for subsequent build_osm_network step
    df_substations.loc[:, "contains"] = df_substations["bus_id"].apply(
        lambda x: x.split("-")[0]
    )

    # Initialise x_node column (if the bus is part of an interconnector) to False, will be set later
    df_substations.loc[:, "x_node"] = False

    # Only included needed columns
    df_substations = df_substations[
        [
            "bus_id",
            "voltage",
            "country",
            "x_node",
            "geometry",
            "polygon",
            "contains",
        ]
    ]

    # Substation data types
    df_substations["voltage"] = df_substations["voltage"].astype(int)

    return df_substations


def _finalise_lines(df_lines):
    """
    Finalises the lines column types.

    Args:
        df_lines (pandas.DataFrame): The input DataFrame containing lines data.

    Returns:
        df_lines (pandas.DataFrame(): The DataFrame with finalised column types
        and transformed data.
    """
    logger.info("Finalising lines column types.")
    df_lines = df_lines.copy()
    # Rename columns
    df_lines.rename(
        columns={
            "id": "line_id",
            "power": "tag_type",
            "frequency": "tag_frequency",
        },
        inplace=True,
    )

    # Initiate new columns for subsequent build_osm_network step
    df_lines.loc[:, "bus0"] = None
    df_lines.loc[:, "bus1"] = None
    df_lines.loc[:, "length"] = None
    df_lines.loc[:, "underground"] = False
    df_lines.loc[df_lines["tag_type"] == "line", "underground"] = False
    df_lines.loc[df_lines["tag_type"] == "cable", "underground"] = True

    # Only include needed columns
    df_lines = df_lines[
        [
            "line_id",
            "circuits",
            "voltage",
            "bus0",
            "bus1",
            "length",
            "underground",
            "geometry",
        ]
    ]

    df_lines["circuits"] = df_lines["circuits"].astype(int)
    df_lines["voltage"] = df_lines["voltage"].astype(int)

    return df_lines


def _finalise_links(df_links):
    """
    Finalises the links column types.

    Args:
        df_links (pandas.DataFrame): The input DataFrame containing links data.

    Returns:
        df_links (pandas.DataFrame(): The DataFrame with finalised column types
        and transformed data.
    """
    logger.info("Finalising links column types.")
    df_links = df_links.copy()
    # Rename columns
    df_links.rename(
        columns={
            "id": "link_id",
            "rating": "p_nom",
        },
        inplace=True,
    )

    # Initiate new columns for subsequent build_osm_network step
    df_links["bus0"] = None
    df_links["bus1"] = None
    df_links["length"] = None
    df_links["underground"] = True

    # Only include needed columns
    df_links = df_links[
        [
            "link_id",
            "voltage",
            "p_nom",
            "bus0",
            "bus1",
            "length",
            "underground",
            "geometry",
        ]
    ]

    df_links["p_nom"] = df_links["p_nom"].astype(int)
    df_links["voltage"] = df_links["voltage"].astype(int)

    return df_links


def _import_substations(path_substations):
    """
    Import substations from the given input paths. This function imports both
    substations from OSM ways as well as relations that contain nested
    information on the substations shape and electrical parameters. Ways and
    relations are subsequently concatenated to form a single DataFrame
    containing unique bus ids.

    Args:
        path_substations (dict): A dictionary containing input paths for
        substations.

    Returns:
        pd.DataFrame: A DataFrame containing the imported substations data.
    """
    cols_substations_way = [
        "id",
        "geometry",
        "country",
        "power",
        "substation",
        "voltage",
        "frequency",
    ]
    cols_substations_relation = [
        "id",
        "country",
        "power",
        "substation",
        "voltage",
        "frequency",
    ]
    df_substations_way = pd.DataFrame(columns=cols_substations_way)
    df_substations_relation = pd.DataFrame(columns=cols_substations_relation)

    logger.info("Importing substations")
    for key in path_substations:
        logger.info(f"Processing {key}...")
        for idx, ip in enumerate(path_substations[key]):
            if (
                os.path.exists(ip) and os.path.getsize(ip) > 400
            ):  # unpopulated OSM json is about 51 bytes
                country = os.path.basename(os.path.dirname(path_substations[key][idx]))
                logger.info(
                    f" - Importing {key} {str(idx+1).zfill(2)}/{str(len(path_substations[key])).zfill(2)}: {ip}"
                )
                with open(ip, "r") as f:
                    data = json.load(f)

                df = pd.DataFrame(data["elements"])
                df["id"] = df["id"].astype(str)
                # new string that adds "way/" to id
                df["id"] = df["id"].apply(
                    lambda x: (
                        f"way/{x}" if key == "substations_way" else f"relation/{x}"
                    )
                )
                df["country"] = country

                col_tags = ["power", "substation", "voltage", "frequency"]

                tags = pd.json_normalize(df["tags"]).map(
                    lambda x: str(x) if pd.notnull(x) else x
                )

                for ct in col_tags:
                    if ct not in tags.columns:
                        tags[ct] = pd.NA

                tags = tags.loc[:, col_tags]

                df = pd.concat([df, tags], axis="columns")

                if key == "substations_way":
                    df.drop(columns=["type", "tags", "bounds", "nodes"], inplace=True)
                    df_substations_way = pd.concat(
                        [df_substations_way, df], axis="rows"
                    )
                elif key == "substations_relation":
                    df.drop(columns=["type", "tags", "bounds"], inplace=True)
                    df_substations_relation = pd.concat(
                        [df_substations_relation, df], axis="rows"
                    )

            else:
                logger.info(
                    f" - Skipping {key} {str(idx+1).zfill(2)}/{str(len(path_substations[key])).zfill(2)} (empty): {ip}"
                )
                continue
        logger.info("---")

    df_substations_way.drop_duplicates(subset="id", keep="first", inplace=True)
    df_substations_relation.drop_duplicates(subset="id", keep="first", inplace=True)

    df_substations_way["geometry"] = df_substations_way.apply(_create_polygon, axis=1)

    # Normalise the members column of df_substations_relation
    cols_members = ["id", "type", "ref", "role", "geometry"]
    df_substations_relation_members = pd.DataFrame(columns=cols_members)

    for index, row in df_substations_relation.iterrows():
        col_members = ["type", "ref", "role", "geometry"]
        df = pd.json_normalize(row["members"])

        for cm in col_members:
            if cm not in df.columns:
                df[cm] = pd.NA

        df = df.loc[:, col_members]
        df["id"] = str(row["id"])
        df["ref"] = df["ref"].astype(str)
        df = df[df["type"] != "node"]
        df = df.dropna(subset=["geometry"])
        df = df[~df["role"].isin(["", "incoming_line", "substation", "inner"])]
        df_substations_relation_members = pd.concat(
            [df_substations_relation_members, df], axis="rows"
        )

    df_substations_relation_members.reset_index(inplace=True)
    df_substations_relation_members["linestring"] = (
        df_substations_relation_members.apply(_create_linestring, axis=1)
    )
    df_substations_relation_members_grouped = (
        df_substations_relation_members.groupby("id")["linestring"]
        .apply(lambda x: linemerge(x.tolist()))
        .reset_index()
    )
    df_substations_relation_members_grouped["geometry"] = (
        df_substations_relation_members_grouped["linestring"].apply(
            lambda x: x.convex_hull
        )
    )

    df_substations_relation = (
        df_substations_relation.join(
            df_substations_relation_members_grouped.set_index("id"), on="id", how="left"
        )
        .drop(columns=["members", "linestring"])
        .dropna(subset=["geometry"])
    )

    # reorder columns and concatenate
    df_substations_relation = df_substations_relation[cols_substations_way]
    df_substations = pd.concat(
        [df_substations_way, df_substations_relation], axis="rows"
    )

    return df_substations


def _remove_lines_within_substations(gdf_lines, gdf_substations_polygon):
    """
    Removes lines that are within substation polygons from the given
    GeoDataFrame of lines. These are not needed to create network (e.g. bus
    bars, switchgear, etc.)

    Parameters:
    - gdf_lines (GeoDataFrame): A GeoDataFrame containing lines with 'line_id'
      and 'geometry' columns.
    - gdf_substations_polygon (GeoDataFrame): A GeoDataFrame containing
      substation polygons.

    Returns:
    GeoDataFrame: A new GeoDataFrame without lines within substation polygons.
    """
    logger.info("Identifying and removing lines within substation polygons...")
    gdf = gpd.sjoin(
        gdf_lines[["line_id", "geometry"]],
        gdf_substations_polygon,
        how="inner",
        predicate="within",
    )["line_id"]

    logger.info(
        f"Removed {len(gdf)} lines within substations of original {len(gdf_lines)} lines."
    )
    gdf_lines = gdf_lines[~gdf_lines["line_id"].isin(gdf)]

    return gdf_lines


def _merge_touching_polygons(df):
    """
    Merge touching polygons in a GeoDataFrame.

    Parameters:
    - df: pandas.DataFrame or geopandas.GeoDataFrame
        The input DataFrame containing the polygons to be merged.

    Returns:
    - gdf: geopandas.GeoDataFrame
        The GeoDataFrame with merged polygons.
    """

    gdf = gpd.GeoDataFrame(df, geometry="polygon", crs=crs)
    combined_polygons = unary_union(gdf.geometry)
    if combined_polygons.geom_type == "MultiPolygon":
        gdf_combined = gpd.GeoDataFrame(
            geometry=[poly for poly in combined_polygons.geoms], crs=crs
        )
    else:
        gdf_combined = gpd.GeoDataFrame(geometry=[combined_polygons], crs=crs)

    gdf.reset_index(drop=True, inplace=True)

    for i, combined_geom in gdf_combined.iterrows():
        mask = gdf.intersects(combined_geom.geometry)
        gdf.loc[mask, "polygon_merged"] = combined_geom.geometry

    gdf.drop(columns=["polygon"], inplace=True)
    gdf.rename(columns={"polygon_merged": "polygon"}, inplace=True)

    return gdf


def _add_endpoints_to_line(linestring, polygon_dict, tol=BUS_TOL / 2):
    """
    Adds endpoints to a line by removing any overlapping areas with polygons.

    Parameters:
    linestring (LineString): The original line to add endpoints to.
    polygon_dict (dict): A dictionary of polygons, where the keys are bus IDs and the values are the corresponding polygons.

    Returns:
    LineString: The modified line with added endpoints.
    """
    if not polygon_dict:
        return linestring
    polygon_pois = {
        bus_id: polylabel(polygon, tol) for bus_id, polygon in polygon_dict.items()
    }
    polygon_unary = unary_union(list(polygon_dict.values()))

    # difference with polygon
    linestring_new = linestring.difference(polygon_unary)

    if type(linestring_new) == MultiLineString:
        # keep the longest line in the multilinestring
        linestring_new = max(linestring_new.geoms, key=lambda x: x.length)

    for p in polygon_pois:
        linestring_new = _add_bus_poi_to_line(linestring_new, polygon_pois[p])

    return linestring_new


def _get_polygons_at_endpoints(linestring, polygon_dict):
    """
    Get the polygons that contain the endpoints of a given linestring.

    Parameters:
    linestring (LineString): The linestring for which to find the polygons at the endpoints.
    polygon_dict (dict): A dictionary containing polygons as values, with bus_ids as keys.

    Returns:
    dict: A dictionary containing bus_ids as keys and polygons as values, where the polygons contain the endpoints of the linestring.
    """
    # Get the endpoints of the linestring
    start_point = Point(linestring.coords[0])
    end_point = Point(linestring.coords[-1])

    # Initialize dictionary to store bus_ids as keys and polygons as values
    bus_id_polygon_dict = {}

    for bus_id, polygon in polygon_dict.items():
        if polygon.contains(start_point) or polygon.contains(end_point):
            bus_id_polygon_dict[bus_id] = polygon

    return bus_id_polygon_dict


def _extend_lines_to_substations(gdf_lines, gdf_substations_polygon, tol=BUS_TOL / 2):
    """
    Extends the lines in the given GeoDataFrame `gdf_lines` to the Pole of Inaccessibility (PoI) of
    the nearest substations represented by the polygons in the
    `gdf_substations_polygon` GeoDataFrame.

    Parameters:
    gdf_lines (GeoDataFrame): A GeoDataFrame containing the lines to be extended.
    gdf_substations_polygon (GeoDataFrame): A GeoDataFrame containing the polygons representing substations.

    Returns:
    GeoDataFrame: A new GeoDataFrame with the lines extended to the substations.
    """
    logger.info(
        "Extending lines ending in substations to interior 'Pole of Inaccessibility'"
    )
    gdf = gpd.sjoin(
        gdf_lines,
        gdf_substations_polygon.drop_duplicates(subset="polygon", inplace=False),
        how="left",
        lsuffix="line",
        rsuffix="bus",
        predicate="intersects",
    ).drop(columns="index_bus")

    # Group by 'line_id' and create a dictionary mapping 'bus_id' to 'geometry_bus', excluding the grouping columns
    gdf = (
        gdf.groupby("line_id")
        .apply(
            lambda x: x[["bus_id", "geometry_bus"]]
            .dropna()
            .set_index("bus_id")["geometry_bus"]
            .to_dict(),
            include_groups=False,
        )
        .reset_index()
    )
    gdf.columns = ["line_id", "bus_dict"]

    gdf["intersects_bus"] = gdf.apply(lambda row: len(row["bus_dict"]) > 0, axis=1)

    gdf.loc[:, "line_geometry"] = gdf.join(
        gdf_lines.set_index("line_id")["geometry"], on="line_id"
    )["geometry"]

    # Polygons at the endpoints of the linestring
    gdf["bus_endpoints"] = gdf.apply(
        lambda row: _get_polygons_at_endpoints(row["line_geometry"], row["bus_dict"]),
        axis=1,
    )

    gdf.loc[:, "line_geometry_new"] = gdf.apply(
        lambda row: _add_endpoints_to_line(
            row["line_geometry"], row["bus_endpoints"], tol
        ),
        axis=1,
    )

    gdf.set_index("line_id", inplace=True)
    gdf_lines.set_index("line_id", inplace=True)

    gdf_lines.loc[:, "geometry"] = gdf["line_geometry_new"]

    return gdf_lines


def _check_if_ways_in_multi(list, longer_list):
    # Check if any of the elements in list are in longer_list
    return any([x in longer_list for x in list])


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_osm_data")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Parameters
    crs = "EPSG:4326"  # Correct crs for OSM data
    min_voltage_ac = 220000  # [unit: V] Minimum voltage value to filter AC lines.
    min_voltage_dc = 150000  #  [unit: V] Minimum voltage value to filter DC links.

    logger.info("---")
    logger.info("SUBSTATIONS")
    # Input
    path_substations = {
        "substations_way": snakemake.input.substations_way,
        "substations_relation": snakemake.input.substations_relation,
    }

    # Cleaning process
    df_substations = _import_substations(path_substations)
    df_substations["voltage"] = _clean_voltage(df_substations["voltage"])

    # Extract converter subset
    df_substations.reset_index(drop=True, inplace=True)
    converter_candidates = (
        df_substations["substation"].str.contains("converter").dropna()
    )
    converter_candidates = converter_candidates[converter_candidates].index
    df_converters = df_substations.loc[converter_candidates]

    df_substations, list_voltages = _filter_by_voltage(
        df_substations, min_voltage=min_voltage_ac
    )
    df_substations["frequency"] = _clean_frequency(df_substations["frequency"])
    df_substations = _clean_substations(df_substations, list_voltages)
    df_substations = _create_substations_geometry(df_substations)

    # Merge touching polygons
    df_substations = _merge_touching_polygons(df_substations)
    df_substations = _create_substations_poi(df_substations)
    df_substations = _finalise_substations(df_substations)

    # Create polygon GeoDataFrame to remove lines within substations
    gdf_substations_polygon = gpd.GeoDataFrame(
        df_substations[["bus_id", "polygon", "voltage"]],
        geometry="polygon",
        crs=crs,
    )

    gdf_substations_polygon["geometry"] = gdf_substations_polygon.polygon.copy()

    # Continue cleaning of converters
    logger.info("---")
    logger.info("CONVERTERS")
    logger.info(f"Extracting {len(df_converters)} converters as subset of substations.")
    df_converters, list_converter_voltages = _filter_by_voltage(
        df_converters, min_voltage=min_voltage_dc
    )
    df_converters.reset_index(drop=True, inplace=True)
    gdf_converters = gpd.GeoDataFrame(
        df_converters[["id", "geometry"]], geometry="geometry", crs=crs
    )

    ### Lines/Cables relations
    logger.info("---")
    logger.info("AC LINES/CABLES RELATIONS")
    path_routes_relation = {
        "routes_relation": snakemake.input.routes_relation,
    }

    df_routes_relation = _import_routes_relation(path_routes_relation)

    df_lines_cables_relation = df_routes_relation.copy()
    df_lines_cables_relation = _drop_duplicate_lines(df_lines_cables_relation)
    df_lines_cables_relation.loc[:, "voltage"] = _clean_voltage(
        df_lines_cables_relation["voltage"]
    )
    df_lines_cables_relation, list_voltages = _filter_by_voltage(
        df_lines_cables_relation, min_voltage=min_voltage_ac
    )
    df_lines_cables_relation.loc[:, "frequency"] = _clean_frequency(
        df_lines_cables_relation["frequency"]
    )
    df_lines_cables_relation = df_lines_cables_relation[
        df_lines_cables_relation["frequency"] != "0"
    ]
    df_lines_cables_relation["frequency"] = "50"
    df_lines_cables_relation.loc[:, "circuits"] = _clean_circuits(
        df_lines_cables_relation["circuits"]
    )
    df_lines_cables_relation.loc[:, "cables"] = _clean_cables(
        df_lines_cables_relation["cables"]
    )
    df_lines_cables_relation = _clean_lines(df_lines_cables_relation, list_voltages)
    df_lines_cables_relation.drop(
        columns=["voltage_original", "cleaned", "circuits_original", "split_elements"],
        inplace=True,
    )

    # Create geometries
    df_lines_cables_relation["components"] = df_lines_cables_relation.apply(
        _create_line, axis=1
    )
    df_lines_cables_relation["geometry"] = df_lines_cables_relation["components"].apply(
        lambda x: x[0]
    )
    df_lines_cables_relation["contains"] = df_lines_cables_relation["components"].apply(
        lambda x: x[1]
    )

    # Show lines that are multilinestring
    subset_multi_linestrings = df_lines_cables_relation[
        df_lines_cables_relation["geometry"].apply(
            lambda x: x.geom_type == "MultiLineString"
        )
    ].copy()
    ways_in_multi_linestring = pd.Series(
        itertools.chain(*subset_multi_linestrings["contains"])
    ).unique()

    # Mark rows that have ways members in subset_multi_linestrings (These cannot be used)
    # TODO: Improve this in the future, to make this constraint less tight (if it is in another multi_linstring), reduce the underlying way/circuit accordingly.
    # TODO: Then, the relation can still be added.
    df_lines_cables_relation["contains_ways_in_multi"] = df_lines_cables_relation[
        "contains"
    ].apply(lambda x: _check_if_ways_in_multi(x, ways_in_multi_linestring))

    # Only keep rows that have contains_ways_in_multi == False
    df_lines_cables_relation = df_lines_cables_relation[
        ~df_lines_cables_relation["contains_ways_in_multi"]
    ]
    # Only keep rows if instance of Linestring
    df_lines_cables_relation = df_lines_cables_relation[
        df_lines_cables_relation["geometry"].apply(
            lambda x: x.geom_type == "LineString"
        )
    ]
    ways_to_replace = pd.Series(
        itertools.chain(*df_lines_cables_relation["contains"])
    ).unique()

    df_lines_cables_relation.rename(columns={"id": "line_id"}, inplace=True)
    df_lines_cables_relation = df_lines_cables_relation[
        ["line_id", "circuits", "voltage", "geometry", "contains"]
    ]
    df_lines_cables_relation["circuits"] = df_lines_cables_relation["circuits"].astype(
        int
    )
    df_lines_cables_relation["voltage"] = df_lines_cables_relation["voltage"].astype(
        int
    )

    # Lines and cables
    logger.info("---")
    logger.info("LINES AND CABLES")
    path_lines = {
        "lines": snakemake.input.lines_way,
        "cables": snakemake.input.cables_way,
    }

    # Import and replace with relations, if relations unique linestrings and line is a member
    df_lines = _import_lines_and_cables(path_lines)
    df_lines = _drop_duplicate_lines(df_lines)

    # Dropping
    len_before = len(df_lines)
    connection_type = (
        df_lines[df_lines["id"].isin(ways_to_replace)].set_index("id")["power"].copy()
    )
    df_lines = df_lines[~df_lines["id"].isin(ways_to_replace)]
    len_after = len(df_lines)
    logger.info(
        f"Dropping {len_before - len_after} OSM ways (AC lines and cables) for their parent OSM relations."
    )

    # Cleaning process
    df_lines.loc[:, "voltage"] = _clean_voltage(df_lines["voltage"])
    df_lines, list_voltages = _filter_by_voltage(df_lines, min_voltage=min_voltage_ac)
    df_lines.loc[:, "circuits"] = _clean_circuits(df_lines["circuits"])
    df_lines.loc[:, "cables"] = _clean_cables(df_lines["cables"])
    df_lines.loc[:, "frequency"] = _clean_frequency(df_lines["frequency"])
    df_lines.loc[:, "wires"] = _clean_wires(df_lines["wires"])

    df_lines = _clean_lines(df_lines, list_voltages)

    # Drop DC lines, will be added through relations later
    len_before = len(df_lines)
    df_lines = df_lines[df_lines["frequency"] == "50"]
    len_after = len(df_lines)
    logger.info(
        f"Dropped {len_before - len_after} DC lines. Keeping {len_after} AC lines."
    )
    df_lines = _create_lines_geometry(df_lines)
    df_lines = _finalise_lines(df_lines)
    df_lines["contains"] = df_lines["line_id"].apply(lambda x: [x.split("-")[0]])

    # Merge
    # Map element in list of column ways to connection_type
    df_lines_cables_relation["underground"] = df_lines_cables_relation[
        "contains"
    ].apply(lambda x: [connection_type[way] for way in x if way in connection_type])
    df_lines_cables_relation["underground"] = df_lines_cables_relation[
        "underground"
    ].apply(lambda x: list(set(x)))
    df_lines_cables_relation["underground"] = df_lines_cables_relation[
        "underground"
    ].apply(lambda x: True if x == ["cable"] else False)

    df_lines = pd.concat(
        [df_lines, df_lines_cables_relation], axis=0, ignore_index=True
    )

    # Create GeoDataFrame
    gdf_lines = gpd.GeoDataFrame(df_lines, geometry="geometry", crs=crs)
    gdf_lines = _remove_lines_within_substations(gdf_lines, gdf_substations_polygon)

    gdf_lines = _extend_lines_to_substations(
        gdf_lines, gdf_substations_polygon, tol=BUS_TOL / 2
    )

    logger.info("---")
    logger.info("HVDC LINKS")

    df_links = df_routes_relation.copy()

    logger.info("Dropping lines without rating.")
    len_before = len(df_links)
    df_links = df_links.dropna(subset=["rating"])
    len_after = len(df_links)
    logger.info(
        f"Dropped {len_before-len_after} elements without rating. "
        + f"Imported {len_after} elements."
    )

    df_links = _drop_duplicate_lines(df_links)
    df_links.loc[:, "voltage"] = _clean_voltage(df_links["voltage"])
    df_links, list_voltages = _filter_by_voltage(df_links, min_voltage=min_voltage_dc)
    # Keep only highest voltage of split string
    df_links.loc[:, "voltage"] = df_links["voltage"].apply(
        lambda x: str(max(map(int, x.split(";"))))
    )
    df_links.loc[:, "frequency"] = _clean_frequency(df_links["frequency"])
    df_links.loc[:, "rating"] = _clean_rating(df_links["rating"])

    df_links.loc[:, "geometry"] = df_links.apply(_create_single_link, axis=1)
    df_links = _finalise_links(df_links)
    gdf_links = gpd.GeoDataFrame(df_links, geometry="geometry", crs=crs).set_index(
        "link_id"
    )

    # Drop polygons and create GDF
    gdf_substations = gpd.GeoDataFrame(
        df_substations.drop(columns=["polygon"]), geometry="geometry", crs=crs
    )

    output_substations_polygon = snakemake.output["substations_polygon"]
    output_substations = snakemake.output["substations"]
    output_converters_polygon = snakemake.output["converters_polygon"]
    output_lines = snakemake.output["lines"]
    output_links = snakemake.output["links"]

    logger.info(
        f"Exporting clean substations with polygon shapes to {output_substations_polygon}"
    )
    gdf_substations_polygon.drop(columns=["geometry"]).to_file(
        output_substations_polygon, driver="GeoJSON"
    )
    logger.info(f"Exporting clean substations to {output_substations}")
    gdf_substations.to_file(output_substations, driver="GeoJSON")
    logger.info(f"Exporting converter polygons to {output_converters_polygon}")
    gdf_converters.to_file(output_converters_polygon, driver="GeoJSON")
    logger.info(f"Exporting clean lines to {output_lines}")
    gdf_lines.to_file(output_lines, driver="GeoJSON")
    logger.info(f"Exporting clean links to {output_links}")
    gdf_links.to_file(output_links, driver="GeoJSON")

    logger.info("Cleaning OSM data completed.")
