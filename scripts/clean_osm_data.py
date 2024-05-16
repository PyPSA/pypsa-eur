# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script is used to clean OpenStreetMap (OSM) data for the PyPSA-Eur 
project.

The script performs various cleaning operations on the OSM data, including:
- Cleaning voltage, circuits, cables, wires, and frequency columns
- Splitting semicolon-separated cells into new rows
- Distributing values to circuits based on the number of splits
- Adding line endings to substations based on line data

The cleaned data is then written to an output file.

Usage:
    python clean_osm_data.py <output_file>

Arguments:
    output_file (str): The path to the output file where the cleaned data will 
    be written.

Example:
    python clean_osm_data.py cleaned_data.csv
"""

import geopandas as gpd
import json
import logging
import os
import numpy as np
import pandas as pd
import re
from shapely.geometry import LineString, Polygon
from shapely.ops import linemerge

from _helpers import configure_logging
logger = logging.getLogger(__name__)


def _create_linestring(row):
    """
    Create a LineString object from the given row.

    Args:
        row (dict): A dictionary containing the row data.

    Returns:
        LineString: A LineString object representing the geometry.

    """
    coords = [(coord['lon'], coord['lat']) for coord in row["geometry"]]
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
    point_coords = [(coord['lon'], coord['lat']) for coord in row["geometry"]]
    
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
        column
        .astype(str)
        .str.lower()
        .str.replace("400/220/110 kV'", "400000;220000;110000")
        .str.replace("400/220/110/20_kv", "400000;220000;110000;20000")
        .str.replace("2x25000", "25000;25000")
    )

    column = (
        column
        .astype(str)
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
    column = column.apply(lambda x: re.sub(r'[^0-9;]', '', str(x)))

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
        column
        .astype(str)
        .str.replace("partial", "")
        .str.replace("1operator=RTE operator:wikidata=Q2178795", "")
        .str.lower()
        .str.replace("1,5", "3")
        .str.replace("1/3", "1")
        .str.replace("<na>", "")
        .str.replace("nan", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r'[^0-9;]', '', x))

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
        column
        .astype(str)
        .str.lower()
        .str.replace("1/3", "1")
        .str.replace("3x2;2", "3")
        .str.replace("<na>", "")
        .str.replace("nan", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r'[^0-9;]', '', x))

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
        column
        .astype(str)
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
    column = column.apply(lambda x: re.sub(r'[^0-9;]', '', x))

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
    voltages = voltage.split(';')
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
        column
        .astype(str)
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
    column = column.apply(lambda x: re.sub(r'[^0-9;.]', '', x))

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
    num_splits = x.groupby('id').size().to_dict()

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


def add_line_endings_tosubstations(substations, lines):
    """
    Add line endings to substations.

    This function takes two pandas DataFrames, `substations` and `lines`, and 
    adds line endings to the substations based on the information from the 
    lines DataFrame.

    Parameters:
    - substations (pandas DataFrame): DataFrame containing information about 
      substations.
    - lines (pandas DataFrame): DataFrame containing information about lines.

    Returns:
    - buses (pandas DataFrame): DataFrame containing the updated information 
      about substations with line endings.

    """
    if lines.empty:
        return substations

    # extract columns from substation df
    bus_s = pd.DataFrame(columns=substations.columns)
    bus_e = pd.DataFrame(columns=substations.columns)

    # Read information from line.csv
    bus_s[["voltage", "country"]] = lines[["voltage", "country"]].astype(str)
    bus_s["geometry"] = lines.geometry.boundary.map(
        lambda p: p.geoms[0] if len(p.geoms) >= 2 else None
    )
    bus_s["lon"] = bus_s["geometry"].map(lambda p: p.x if p != None else None)
    bus_s["lat"] = bus_s["geometry"].map(lambda p: p.y if p != None else None)
    bus_s["bus_id"] = (
        (substations["bus_id"].max() if "bus_id" in substations else 0)
        + 1
        + bus_s.index
    )
    bus_s["dc"] = lines["dc"]

    bus_e[["voltage", "country"]] = lines[["voltage", "country"]].astype(str)
    bus_e["geometry"] = lines.geometry.boundary.map(
        lambda p: p.geoms[1] if len(p.geoms) >= 2 else None
    )
    bus_e["lon"] = bus_e["geometry"].map(lambda p: p.x if p != None else None)
    bus_e["lat"] = bus_e["geometry"].map(lambda p: p.y if p != None else None)
    bus_e["bus_id"] = bus_s["bus_id"].max() + 1 + bus_e.index
    bus_e["dc"] = lines["dc"]

    bus_all = pd.concat([bus_s, bus_e], ignore_index=True)

    # Initialize default values
    bus_all["station_id"] = np.nan
    # Assuming substations completed for installed lines
    bus_all["under_construction"] = False
    bus_all["tag_area"] = 0.0
    bus_all["symbol"] = "substation"
    # TODO: this tag may be improved, maybe depending on voltage levels
    bus_all["tag_substation"] = "transmission"
    bus_all["tag_source"] = "line_ending"

    buses = pd.concat([substations, bus_all], ignore_index=True)

    # # Assign index to bus_id
    buses["bus_id"] = buses.index

    # TODO: pypsa-eur: change this later to improve country assignment
    bool_multiple_countries = buses["country"].str.contains(";")
    buses.loc[bool_multiple_countries, "country"] = buses.loc[bool_multiple_countries, "country"].str.split(";").str[0]

    return buses


def _import_lines_and_cables(input_path_lines_cables):
    """
    Import lines and cables from the given input paths.

    Parameters:
    - input_path_lines_cables (dict): A dictionary containing the input paths for lines and cables data.

    Returns:
    - df_lines (DataFrame): A DataFrame containing the imported lines and cables data.

    """
    columns = ["id", "bounds", "nodes", "geometry", "country", "power", "cables", "circuits", "frequency", "voltage", "wires"]
    df_lines = pd.DataFrame(columns=columns)

    logger.info("Importing lines and cables")
    for key in input_path_lines_cables:
        logger.info(f"Processing {key}...")
        for idx, ip in enumerate(input_path_lines_cables[key]):
            if os.path.exists(ip) and os.path.getsize(ip) > 400: # unpopulated OSM json is about 51 bytes
                country = os.path.basename(os.path.dirname(input_path_lines_cables[key][idx]))
                
                logger.info(f" - Importing {key} {str(idx+1).zfill(2)}/{str(len(input_path_lines_cables[key])).zfill(2)}: {ip}")
                with open(ip, "r") as f:
                    data = json.load(f)
                
                df = pd.DataFrame(data['elements'])
                df["id"] = df["id"].astype(str)
                df["country"] = country

                col_tags = ["power", "cables", "circuits", "frequency", "voltage", "wires"]

                tags = pd.json_normalize(df["tags"]) \
                    .map(lambda x: str(x) if pd.notnull(x) else x)
                
                for ct in col_tags:
                    if ct not in tags.columns:
                        tags[ct] = pd.NA
                
                tags = tags.loc[:, col_tags]

                df = pd.concat([df, tags], axis="columns") 
                df.drop(columns=["type", "tags"], inplace=True)
                
                df_lines = pd.concat([df_lines, df], axis="rows")

            else:
                logger.info(f" - Skipping {key} {str(idx+1).zfill(2)}/{str(len(input_path_lines_cables[key])).zfill(2)} (empty): {ip}")
                continue
        logger.info("---")
    
    return df_lines


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
    duplicate_rows = df_lines[df_lines.duplicated(subset=['id'], keep=False)].copy()

    # Group rows by id and aggregate the country column to a string split by semicolon
    grouped_duplicates = duplicate_rows.groupby('id')["country"].agg(lambda x: ';'.join(x)).reset_index()
    duplicate_rows.drop_duplicates(subset="id", inplace=True)
    duplicate_rows.drop(columns=["country"], inplace=True)
    duplicate_rows = duplicate_rows.join(grouped_duplicates.set_index('id'), on='id', how='left')

    # Drop duplicates and update the df_lines dataframe with the cleaned data
    df_lines = df_lines[~df_lines["id"].isin(duplicate_rows["id"])]
    df_lines = pd.concat([df_lines, duplicate_rows], axis="rows")

    return df_lines


def _filter_lines_by_voltage(df_lines, voltage_min=200000):
    """
    Filter lines in the DataFrame `df_lines` based on the voltage in V.

    Parameters:
    - df_lines (pandas.DataFrame): The DataFrame containing the lines data.
    - voltage_min (int, optional): The minimum voltage value to filter the 
      lines. Defaults to 200000 [unit: V].

    Returns:
    - filtered df_lines (pandas.DataFrame): The filtered DataFrame containing 
      the lines data above voltage_min.
    - list_voltages (list): A list of unique voltage values above voltage_min.
      The type of the list elements is string.
    """
    logger.info(f"Filtering lines by voltage. Only keeping lines above and including {voltage_min} V.")
    list_voltages = df_lines["voltage"].str.split(";").explode().unique().astype(str)
    # Keep numeric strings
    list_voltages = list_voltages[np.vectorize(str.isnumeric)(list_voltages)]
    list_voltages = list_voltages.astype(int)
    list_voltages = list_voltages[list_voltages >= int(voltage_min)]
    list_voltages = list_voltages.astype(str)

    bool_voltages = df_lines["voltage"].apply(_check_voltage, list_voltages=list_voltages)
    df_lines = df_lines[bool_voltages]

    return df_lines, list_voltages


def _clean_lines(df_lines):
    """
    Cleans and processes the `df_lines` DataFrame heuristically based on the 
    information available per respective line and cable.
    Further checks to ensure data consistency and completeness.

    Parameters
    ----------
    df_lines : pandas.DataFrame
        The input DataFrame containing line information with columns such as 
        'voltage', 'circuits', 'frequency', 'cables', 'split_elements', 'id', 
        etc.

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
    bool_voltages = df_lines["voltage"].apply(_check_voltage, list_voltages=list_voltages)
    df_lines = df_lines[bool_voltages]

    bool_ac = df_lines["frequency"] != "0"
    bool_dc = ~bool_ac
    valid_frequency = ["50", "0"]
    bool_invalid_frequency = df_lines["frequency"].apply(lambda x: x not in valid_frequency)

    bool_noinfo = (df_lines["cables"] == "") & (df_lines["circuits"] == "")
    # Fill in all values where cables info and circuits does not exist. Assuming 1 circuit
    df_lines.loc[bool_noinfo, "circuits"] = "1"
    df_lines.loc[bool_noinfo & bool_invalid_frequency, "frequency"] = "50"
    df_lines.loc[bool_noinfo, "cleaned"] = True

    # Fill in all values where cables info exists and split_elements == 1
    bool_cables_ac = (df_lines["cables"] != "") & \
        (df_lines["split_elements"] == 1) & \
        (df_lines["cables"] != "0") & \
        (df_lines["cables"].apply(lambda x: len(x.split(";")) == 1)) & \
        (df_lines["circuits"] == "") & \
        (df_lines["cleaned"] == False) & \
        bool_ac
    
    df_lines.loc[bool_cables_ac, "circuits"] = df_lines.loc[bool_cables_ac, "cables"] \
        .apply(lambda x: str(int(max(1, np.floor_divide(int(x),3)))))
    
    df_lines.loc[bool_cables_ac, "frequency"] = "50"
    df_lines.loc[bool_cables_ac, "cleaned"] = True

    bool_cables_dc = (df_lines["cables"] != "") & \
        (df_lines["split_elements"] == 1) & \
        (df_lines["cables"] != "0") & \
        (df_lines["cables"].apply(lambda x: len(x.split(";")) == 1)) & \
        (df_lines["circuits"] == "") & \
        (df_lines["cleaned"] == False) & \
        bool_dc
    
    df_lines.loc[bool_cables_dc, "circuits"] = df_lines.loc[bool_cables_dc, "cables"] \
        .apply(lambda x: str(int(max(1, np.floor_divide(int(x),2)))))
    
    df_lines.loc[bool_cables_dc, "frequency"] = "0"
    df_lines.loc[bool_cables_dc, "cleaned"] = True

    # Fill in all values where circuits info exists and split_elements == 1
    bool_lines = (df_lines["circuits"] != "") & \
        (df_lines["split_elements"] == 1) & \
        (df_lines["circuits"] != "0") & \
        (df_lines["circuits"].apply(lambda x: len(x.split(";")) == 1)) & \
        (df_lines["cleaned"] == False) 
    
    df_lines.loc[bool_lines & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_lines & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_lines, "cleaned"] = True

    # Clean those values where number of voltages split by semicolon is larger than no cables or no circuits
    bool_cables = (df_lines["voltage_original"].apply(lambda x: len(x.split(";")) > 1)) & \
        (df_lines["cables"].apply(lambda x: len(x.split(";")) == 1)) & \
        (df_lines["circuits"].apply(lambda x: len(x.split(";")) == 1)) & \
        (df_lines["cleaned"] == False)
    
    df_lines.loc[bool_cables, "circuits"] = df_lines[bool_cables] \
        .apply(_distribute_to_circuits, axis=1)
    df_lines.loc[bool_cables & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_cables & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_cables, "cleaned"] = True

    # Clean those values where multiple circuit values are present, divided by semicolon
    bool_cables = (df_lines["circuits"].apply(lambda x: len(x.split(";")) > 1)) & \
        (df_lines.apply(lambda row: len(row["circuits"].split(";")) == row["split_elements"], axis=1)) & \
        (df_lines["cleaned"] == False)
    
    df_lines.loc[bool_cables, "circuits"] = df_lines.loc[bool_cables] \
        .apply(lambda row: str(row["circuits"].split(";")[
            int(row["id"].split("-")[-1])-1
        ]), axis=1)
    
    df_lines.loc[bool_cables & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_cables & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_cables, "cleaned"] = True

    # Clean those values where multiple cables values are present, divided by semicolon
    bool_cables = (df_lines["cables"].apply(lambda x: len(x.split(";")) > 1)) & \
        (df_lines.apply(lambda row: len(row["cables"].split(";")) == row["split_elements"], axis=1)) & \
        (df_lines["cleaned"] == False)

    df_lines.loc[bool_cables, "circuits"] = df_lines.loc[bool_cables] \
        .apply(lambda row: 
            str(max(1,
                np.floor_divide(
                    int(row["cables"].split(";")[int(row["id"].split("-")[-1])-1]),
                    3
                    )
                )),
            axis=1)
    
    df_lines.loc[bool_cables & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_cables & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_cables, "cleaned"] = True

    # All remaining lines to circuits == 1
    bool_leftover = (df_lines["cleaned"] == False)
    if sum(bool_leftover) > 0:
        str_id = "; ".join(str(id) for id in df_lines.loc[bool_leftover, "id"])
        logger.info(f"Setting circuits of remaining {sum(bool_leftover)} lines to 1...")
        logger.info(f"Lines affected: {str_id}")
    
    df_lines.loc[bool_leftover, "circuits"] = "1"
    df_lines.loc[bool_leftover & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_leftover & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_leftover, "cleaned"] = True

    return df_lines


def _finalise_lines(df_lines):
    """
    Finalises the lines column types and creates geometries.

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
            "frequency":"tag_frequency",
            }, inplace=True)
    
    # Initiate new columns for subsequent build_osm_network step
    df_lines.loc[:, "bus0"] = None
    df_lines.loc[:, "bus1"] = None
    df_lines.loc[:, "length"] = None
    df_lines.loc[:, "underground"] = False
    df_lines.loc[df_lines["tag_type"] == "line", "underground"] = False
    df_lines.loc[df_lines["tag_type"] == "cable", "underground"] = True
    df_lines.loc[:, "under_construction"] = False
    df_lines.loc[:, "dc"] = False
    df_lines.loc[df_lines["tag_frequency"] == "50", "dc"] = False
    df_lines.loc[df_lines["tag_frequency"] == "0", "dc"] = True

    # Only include needed columns
    df_lines = df_lines[[
        "line_id",
        "circuits",
        "tag_type",
        "voltage",
        "tag_frequency",
        "bus0",
        "bus1",
        "length",
        "underground",
        "under_construction",
        "dc",
        "country",
        "geometry",
        ]]
    
    # Set lines data types df.apply(pd.to_numeric, args=('coerce',))
    # This workaround is needed as otherwise the column dtypes remain "objects"
    df_lines.loc[:, "circuits_num"] = df_lines["circuits"].astype(int)
    df_lines.loc[:, "voltage_num"] = df_lines["voltage"].astype(int)
    df_lines.loc[:, "tag_frequency_num"] = df_lines["tag_frequency"].astype(int)
    df_lines.drop(columns=["circuits", "voltage", "tag_frequency"], inplace=True)

    col_rename_dict = {
        "circuits_num": "circuits",
        "voltage_num": "voltage",
        "tag_frequency_num": "tag_frequency"
    }   

    df_lines.rename(columns=col_rename_dict, inplace=True)

    # Create shapely linestrings from geometries
    df_lines.loc[:, "geometry"] = df_lines.apply(_create_linestring, axis=1)  

    # Drop all rows where the geometry has equal start and end point
    # These are usually not lines, but outlines of areas.
    bool_circle = df_lines["geometry"].apply(lambda x: x.coords[0] == x.coords[-1]) 
    df_lines = df_lines[~bool_circle] 

    return df_lines


def _import_substations(input_path_substations):
    """
    Import substations from the given input paths. This function imports both
    substations from OSM ways as well as relations that contain nested 
    information on the substations shape and electrical parameters. Ways and
    relations are subsequently concatenated to form a single DataFrame 
    containing unique bus ids.

    Args:
        input_path_substations (dict): A dictionary containing input paths for 
        substations.

    Returns:
        pd.DataFrame: A DataFrame containing the imported substations data.
    """
    cols_substations_way = ["id", "geometry", "country", "power", "substation", "voltage", "frequency"]
    cols_substations_relation = ["id", "country", "power", "substation", "voltage", "frequency"]
    df_substations_way = pd.DataFrame(columns = cols_substations_way)
    df_substations_relation = pd.DataFrame(columns = cols_substations_relation)

    logger.info("Importing substations")
    for key in input_path_substations:
        logger.info(f"Processing {key}...")
        for idx, ip in enumerate(input_path_substations[key]):
            if os.path.exists(ip) and os.path.getsize(ip) > 400: # unpopulated OSM json is about 51 bytes
                country = os.path.basename(os.path.dirname(input_path_substations[key][idx]))  
                logger.info(f" - Importing {key} {str(idx+1).zfill(2)}/{str(len(input_path_substations[key])).zfill(2)}: {ip}")
                with open(ip, "r") as f:
                    data = json.load(f)
                
                df = pd.DataFrame(data['elements'])
                df["id"] = df["id"].astype(str)
                # new string that adds "way/" to id
                df["id"] = df["id"].apply(lambda x: f"way/{x}" if key == "substations_way" else f"relation/{x}")
                df["country"] = country

                col_tags = ["power", "substation", "voltage", "frequency"]

                tags = pd.json_normalize(df["tags"]) \
                    .map(lambda x: str(x) if pd.notnull(x) else x)
                
                for ct in col_tags:
                    if ct not in tags.columns:
                        tags[ct] = pd.NA
                
                tags = tags.loc[:, col_tags]

                df = pd.concat([df, tags], axis="columns") 

                if key == "substations_way":
                    df.drop(columns=["type", "tags", "bounds", "nodes"], inplace=True)
                    df_substations_way = pd.concat([df_substations_way, df], axis="rows")
                elif key == "substations_relation":
                    df.drop(columns=["type", "tags", "bounds"], inplace=True)
                    df_substations_relation = pd.concat([df_substations_relation, df], axis="rows")

            else:
                logger.info(f" - Skipping {key} {str(idx+1).zfill(2)}/{str(len(input_path_substations[key])).zfill(2)} (empty): {ip}")
                continue
        logger.info("---")

    df_substations_way.drop_duplicates(subset='id', keep='first', inplace=True)
    df_substations_relation.drop_duplicates(subset='id', keep='first', inplace=True)

    df_substations_way["geometry"] = df_substations_way.apply(_create_polygon, axis=1)

    # Normalise the members column of df_substations_relation
    cols_members = ["id", "type", "ref", "role", "geometry"]
    df_substations_relation_members = pd.DataFrame(columns = cols_members)

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
        df_substations_relation_members = pd.concat([df_substations_relation_members, df], axis="rows")
    
    df_substations_relation_members.reset_index(inplace=True)
    df_substations_relation_members["linestring"] = df_substations_relation_members.apply(_create_linestring, axis=1)  
    df_substations_relation_members_grouped = df_substations_relation_members.groupby('id')['linestring'] \
        .apply(lambda x: linemerge(x.tolist())).reset_index()
    df_substations_relation_members_grouped["geometry"] = df_substations_relation_members_grouped["linestring"].apply(lambda x: x.convex_hull)
    
    df_substations_relation = df_substations_relation.join(
        df_substations_relation_members_grouped.set_index('id'), 
        on='id', how='left'
        ).drop(columns=["members", "linestring"]) \
        .dropna(subset=["geometry"])
    
    # reorder columns and concatenate
    df_substations_relation = df_substations_relation[cols_substations_way]
    df_substations = pd.concat([df_substations_way, df_substations_relation], axis="rows")

    return df_substations

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_osm_data")
    
    configure_logging(snakemake)
    
    # Parameters
    crs = "EPSG:4326"       # Correct crs for OSM data
    voltage_min = 200000    # [unit: V] Minimum voltage value to filter lines. 
    
    # TODO pypsa-eur: Temporary solution as one AC line between converters will 
    # create an error in simplify_network:
    lines_to_drop = ["775580659"]

    # Input
    input_path_substations = {
        "substations_way": snakemake.input.substations_way,
        "substations_relation": snakemake.input.substations_relation,
    }

    input_path_lines_cables = {
        "lines": snakemake.input.lines_way,
        "cables": snakemake.input.cables_way,
    }

    # Cleaning process
    df_lines = _import_lines_and_cables(input_path_lines_cables)
    df_lines = _drop_duplicate_lines(df_lines)
    df_lines.loc[:, "voltage"] = _clean_voltage(df_lines["voltage"])
    df_lines, list_voltages = _filter_lines_by_voltage(df_lines, voltage_min=voltage_min)

    df_lines.loc[:, "circuits"] = _clean_circuits(df_lines["circuits"])
    df_lines.loc[:, "cables"] = _clean_cables(df_lines["cables"])
    df_lines.loc[:, "frequency"] = _clean_frequency(df_lines["frequency"])
    df_lines.loc[:, "wires"] = _clean_wires(df_lines["wires"])

    df_lines = _clean_lines(df_lines)
    df_lines = _finalise_lines(df_lines)
    
    # Dropping specific lines, manually
    if lines_to_drop in df_lines["line_id"].values:
        df_lines.drop(df_lines[df_lines["line_id"].isin(lines_to_drop)].index, inplace=True)
    
    # Create GeoDataFrame
    gdf_lines = gpd.GeoDataFrame(df_lines, geometry = "geometry", crs = crs)

    ############# BUSES / SUBSTATIONS ######################
    df_substations = _import_substations(input_path_substations)
 

    # Create centroids from geometries
    df_substations.loc[:, "polygon"] = df_substations["geometry"]
    df_substations.loc[:, "geometry"] = df_substations["geometry"].apply(lambda x: x.centroid)
    df_substations.loc[:, "lon"] = df_substations["geometry"].apply(lambda x: x.x)
    df_substations.loc[:, "lat"] = df_substations["geometry"].apply(lambda x: x.y)

    logger.info("Cleaning substations")
    # Clean columns
    df_substations["voltage"] = _clean_voltage(df_substations["voltage"])
    df_substations["frequency"] = _clean_frequency(df_substations["frequency"])
    df_substations["frequency"] = df_substations["frequency"].astype(str, errors="ignore")

    list_voltages = df_substations["voltage"].str.split(";").explode().unique().astype(str)
    list_voltages = list_voltages[np.vectorize(len)(list_voltages) >= 6]
    list_voltages = list_voltages[~np.char.startswith(list_voltages, '1')]

    bool_voltages = df_substations["voltage"].apply(_check_voltage, list_voltages=list_voltages)
    df_substations = df_substations[bool_voltages]

    df_substations = _split_cells(df_substations)
    bool_voltages = df_substations["voltage"].apply(_check_voltage, list_voltages=list_voltages)
    df_substations = df_substations[bool_voltages]
    df_substations["split_count"] = df_substations["id"].apply(lambda x: x.split("-")[1] if "-" in x else "0")
    df_substations["split_count"] = df_substations["split_count"].astype(int)

    bool_split = df_substations["split_elements"] > 1
    bool_frequency_len = df_substations["frequency"].apply(lambda x: len(x.split(";"))) == df_substations["split_elements"]
    df_substations.loc[bool_frequency_len & bool_split, "frequency"] = df_substations.loc[bool_frequency_len & bool_split, "frequency"] \
    
    op_freq = lambda row: row["frequency"].split(";")[row["split_count"]-1]

    df_substations.loc[bool_frequency_len & bool_split, ["frequency"]] = df_substations.loc[bool_frequency_len & bool_split, ] \
        .apply(op_freq, axis=1)
    
    df_substations = _split_cells(df_substations, cols=["frequency"])
    bool_invalid_frequency = df_substations["frequency"].apply(lambda x: x not in ["50", "0"])
    df_substations.loc[bool_invalid_frequency, "frequency"] = "50"
    df_substations["power"] = "substation"
    df_substations["substation"] = "transmission"
    df_substations["dc"] = False
    df_substations.loc[df_substations["frequency"] == "0", "dc"] = True
    df_substations["under_construction"] = False
    df_substations["station_id"] = None
    df_substations["tag_area"] = None
    df_substations["tag_source"] = df_substations["id"]

    gdf_substations_polygon = gpd.GeoDataFrame(
        df_substations[["id", "polygon"]], 
        geometry = "polygon", 
        crs = "EPSG:4326"
        )
    
    filepath_substations_polygon = snakemake.output["substations_polygon"]
    # save substations output
    logger.info(f"Exporting clean substations with polygon shapes to {filepath_substations_polygon}")
    parentfolder_substations_polygon = os.path.dirname(filepath_substations_polygon)
    if not os.path.exists(parentfolder_substations_polygon):
        # Create the folder and its parent directories if they don't exist
        os.makedirs(parentfolder_substations_polygon)

    logger.info(f"Exporting clean substations to {filepath_substations_polygon}")
    gdf_substations_polygon.to_file(filepath_substations_polygon, driver="GeoJSON")    
    

    logger.info("Identifying and removing lines within substation polygons...")
    lines_within_substations = gpd.sjoin(
        gdf_lines[["line_id", "geometry"]], 
        gdf_substations_polygon, 
        how = "inner",
        predicate = "within"
        )["line_id"]

    logger.info(f"Removed {len(lines_within_substations)}/{len(gdf_lines)} lines within substations.")
    gdf_lines = gdf_lines[~gdf_lines["line_id"].isin(lines_within_substations)]
    
    # # Create an empty list to store the results
    # results = []

    # subset a to find only country equal to "BE"
    # a[a["country"] == "BE"]

    # logger.info("Identifying and removing lines within substation polygons...")
    # for index, row in tqdm(gdf_lines.iterrows(), total=len(gdf_lines)):
    #     line = row['geometry']  
    #     # Check if the LineString is within any Polygon in 'substations_df'
    #     is_within_any_substation = any(line.within(substation_polygon) for substation_polygon in df_substations["polygon"])
    #     results.append(is_within_any_substation)

    # # Add the results to 'gdf_lines'
    # gdf_lines['within_substation'] = results

    # gdf_lines = gdf_lines[~gdf_lines["within_substation"]]
    # logger.info(f"Removed {sum(results)} lines within substations.")

    filepath_lines = snakemake.output["lines"]
    # save substations output
    logger.info(f"Exporting clean lines to {filepath_lines}")
    parentfolder_lines = os.path.dirname(filepath_lines)
    if not os.path.exists(parentfolder_lines):
        # Create the folder and its parent directories if they don't exist
        os.makedirs(parentfolder_lines)

    logger.info(f"Exporting clean lines to {filepath_lines}")
    gdf_lines.to_file(filepath_lines, driver="GeoJSON")

    # rename columns
    df_substations.rename(
        columns={
            "id": "bus_id", 
            "power": "symbol",
            "substation":"tag_substation",
            }, inplace=True)
    
    df_substations = df_substations[[
        "bus_id",
        "symbol", 
        "tag_substation", 
        "voltage", 
        "lon", 
        "lat", 
        "dc", 
        "under_construction", 
        "station_id", 
        "tag_area", 
        "country",
        "geometry",
        "tag_source",
        ]]
    
    df_substations["bus_id"] = df_substations.index

    logger.info("Adding line endings to substations")
    df_substations = add_line_endings_tosubstations(
                df_substations, gdf_lines
            )
    
    #group gdf_substations by voltage and and geometry (dropping duplicates)
    df_substations = df_substations.groupby(["voltage", "lon", "lat", "dc", "tag_source"]).first().reset_index()
    df_substations["bus_id"] = df_substations.index
    
    gdf_substations = gpd.GeoDataFrame(df_substations, geometry = "geometry", crs = "EPSG:4326")

    # Substation data types
    gdf_substations["bus_id"] = gdf_substations["bus_id"].astype(int)
    gdf_substations["voltage"] = gdf_substations["voltage"].astype(int)

    filepath_substations = snakemake.output["substations"]
    # save substations output
    logger.info(f"Exporting clean substations to {filepath_substations}")
    parentfolder_substations = os.path.dirname(filepath_substations)
    if not os.path.exists(parentfolder_substations):
        # Create the folder and its parent directories if they don't exist
        os.makedirs(parentfolder_substations)

    logger.info(f"Exporting clean substations to {filepath_substations}")
    gdf_substations.to_file(filepath_substations, driver="GeoJSON")    
    