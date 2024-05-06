# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
TODO To fill later
"""

from branca.element import Figure
import folium
import geopandas as gpd
import json
import logging
import os
import numpy as np
import pandas as pd
import re
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import linemerge
import tqdm.auto as tqdm

from _helpers import configure_logging
logger = logging.getLogger(__name__)

def clean_osm_data(output):
    with open(output, "w") as file:
        file.write("Hello, world!\n")


def _create_linestring(row):
    coords = [(coord['lon'], coord['lat']) for coord in row["geometry"]]
    return LineString(coords)


def _create_polygon(row):
    """
    Create a Shapely Polygon from a list of coordinate dictionaries.
    
    Parameters:
        coords (list): List of dictionaries with 'lat' and 'lon' keys representing coordinates.
        
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
    Function to clean the raw circuits column: manual fixing and drop nan values

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
    column = column.copy()
    column = (
        column
        .astype(str)
        .str.replace("partial", "")
        .str.replace("1operator=RTE operator:wikidata=Q2178795", "")
        .str.lower()
        .str.replace("1,5", "3") # (way 998005838, should be corrected in OSM soon)
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


def _set_frequency(column):
    column = column.copy()
    to_fifty = column.astype(str) != "0"
    column[to_fifty] = "50"    

    return column


def _check_voltage(voltage, list_voltages):
    voltages = voltage.split(';')
    for v in voltages:
        if v in list_voltages:
            return True
    return False


def _clean_frequency(column):   
    column = column.copy()
    """
    Function to clean the raw frequency column: manual fixing and drop nan values

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
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
    # Create a dictionary to store the number of splits associated with each original ID
    num_splits = {}

    # Split cells and create new rows
    x = df.assign(**{col: df[col].str.split(";") for col in cols})
    x = x.explode(cols, ignore_index=True)

    # Count the number of splits associated with each original ID
    num_splits = x.groupby('id').size().to_dict()

    # Update the 'split_elements' column
    x["split_elements"] = x["id"].map(num_splits)

    # Function to generate the new ID with suffix and update the number of splits
    def generate_new_id(row):
        original_id = row["id"]
        if row["split_elements"] == 1:
            return original_id
        else:
            suffix_counts[original_id] = suffix_counts.get(original_id, 0) + 1
            return f"{original_id}_{suffix_counts[original_id]}"

    # Update the ID column with the new IDs
    x["id"] = x.apply(generate_new_id, axis=1)

    return x


def _distribute_to_circuits(row):
    if row["circuits"] != "":
        circuits = int(row["circuits"])
    else:
        cables = int(row["cables"])
        circuits = cables / 3

    single_circuit = int(max(1, np.floor_divide(circuits, row["split_elements"])))
    single_circuit = str(single_circuit)

    return single_circuit


# Function to check if any substring is in valid_strings
def _any_substring_in_list(s, list_strings):
    substrings = s.split(';')
    return any(sub in list_strings for sub in substrings)


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_osm_data")
    
    configure_logging(snakemake)
    logger.info("Dummy log: clean_osm_data()")

    ############# BUSES / SUBSTATIONS ######################
    input_path_substations = {
        "substations_way": snakemake.input.substations_way,
        "substations_relation": snakemake.input.substations_relation,
    }

    cols_substations_way = ["id", "geometry", "country", "power", "substation", "voltage", "frequency"]
    cols_substations_relation = ["id", "country", "power", "substation", "voltage", "frequency"]
    df_substations_way = pd.DataFrame(columns = cols_substations_way)
    df_substations_relation = pd.DataFrame(columns = cols_substations_relation)

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

    # Create centroids from geometries
    df_substations.loc[:, "geometry"] = df_substations["geometry"].apply(lambda x: x.centroid)
    df_substations.loc[:, "lon"] = df_substations["geometry"].apply(lambda x: x.x)
    df_substations.loc[:, "lat"] = df_substations["geometry"].apply(lambda x: x.y)

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
    df_substations["split_count"] = df_substations["id"].apply(lambda x: x.split("_")[1] if "_" in x else "0")
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
        ]]
    
    gdf_substations = gpd.GeoDataFrame(df_substations, geometry = "geometry", crs = "EPSG:4326")

    filepath_substations = snakemake.output["substations"]
    # save substations output
    logger.info(f"Exporting clean substations to {filepath_substations}")
    parentfolder_substations = os.path.dirname(filepath_substations)
    if not os.path.exists(parentfolder_substations):
        # Create the folder and its parent directories if they don't exist
        os.makedirs(parentfolder_substations)

    gdf_substations.to_file(filepath_substations, driver="GeoJSON")

    ############# LINES AND CABLES ######################

    input_path_lines_cables = {
        "lines": snakemake.input.lines_way,
        "cables": snakemake.input.cables_way,
    }

    columns = ["id", "bounds", "nodes", "geometry", "country", "power", "cables", "circuits", "frequency", "voltage", "wires"]
    df_lines = pd.DataFrame(columns=columns)
    crs = "EPSG:4326"

    # using tqdm loop over input path

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

    # Find duplicates based on id column
    duplicate_rows = df_lines[df_lines.duplicated(subset=['id'], keep=False)].copy()
    # group rows by id and aggregate the country column to a string split by semicolon
    grouped_duplicates = duplicate_rows.groupby('id')["country"].agg(lambda x: ';'.join(x)).reset_index()
    duplicate_rows.drop_duplicates(subset="id", inplace=True)
    duplicate_rows.drop(columns=["country"], inplace=True)
    duplicate_rows = duplicate_rows.join(grouped_duplicates.set_index('id'), on='id', how='left')

    # Drop duplicates and update the df_lines dataframe with the cleaned data
    df_lines = df_lines[~df_lines["id"].isin(duplicate_rows["id"])]
    df_lines = pd.concat([df_lines, duplicate_rows], axis="rows")

    # Initiate boolean with False, only set to true if all cleaning steps are passed
    df_lines["cleaned"] = False
    df_lines["voltage"] = _clean_voltage(df_lines["voltage"])

    list_voltages = df_lines["voltage"].str.split(";").explode().unique().astype(str)
    list_voltages = list_voltages[np.vectorize(len)(list_voltages) >= 6]
    list_voltages = list_voltages[~np.char.startswith(list_voltages, '1')]

    bool_voltages = df_lines["voltage"].apply(_check_voltage, list_voltages=list_voltages)
    df_lines = df_lines[bool_voltages]

    # Additional cleaning
    df_lines["circuits"] = _clean_circuits(df_lines["circuits"])
    df_lines["cables"] = _clean_cables(df_lines["cables"])
    df_lines["frequency"] = _clean_frequency(df_lines["frequency"])
    df_lines["wires"] = _clean_wires(df_lines["wires"])

    df_lines["voltage_original"] = df_lines["voltage"]
    df_lines["circuits_original"] = df_lines["circuits"]

    df_lines = _split_cells(df_lines)
    bool_voltages = df_lines["voltage"].apply(_check_voltage, list_voltages=list_voltages)
    df_lines = df_lines[bool_voltages]

    bool_ac = df_lines["frequency"] != "0"
    bool_dc = ~bool_ac
    bool_noinfo = (df_lines["cables"] == "") & (df_lines["circuits"] == "")
    valid_frequency = ["50", "0"]
    bool_invalid_frequency = df_lines["frequency"].apply(lambda x: x not in valid_frequency)

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
            int(row["id"].split("_")[-1])-1
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
                    int(row["cables"].split(";")[int(row["id"].split("_")[-1])-1]),
                    3
                    )
                )),
            axis=1)
    
    df_lines.loc[bool_cables & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_cables & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_cables, "cleaned"] = True

    # All remaining lines to circuits == 1
    bool_leftover = (df_lines["cleaned"] == False)
    str_id = "; ".join(str(id) for id in df_lines.loc[bool_leftover, "id"])
    logger.info(f"Setting circuits of remaining {sum(bool_leftover)} lines to 1...")
    logger.info(f"Lines affected: {str_id}")
    df_lines.loc[bool_leftover, "circuits"] = "1"
    df_lines.loc[bool_leftover & bool_ac, "frequency"] = "50"
    df_lines.loc[bool_leftover & bool_dc, "frequency"] = "0"
    df_lines.loc[bool_leftover, "cleaned"] = True

    # rename columns
    df_lines.rename(
        columns={
            "id": "line_id", 
            "power": "tag_type",
            "frequency":"tag_frequency",
            }, inplace=True)
    
    df_lines["bus0"] = None
    df_lines["bus1"] = None
    df_lines["length"] = None
    df_lines.loc[df_lines["tag_type"] == "line", "underground"] = False
    df_lines.loc[df_lines["tag_type"] == "cable", "underground"] = True
    df_lines["under_construction"] = False
    df_lines.loc[df_lines["tag_frequency"] == "0", "dc"] = True
    df_lines.loc[df_lines["tag_frequency"] == "50", "dc"] = False

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
    
    df_lines["geometry"] = df_lines.apply(_create_linestring, axis=1)  
    gdf_lines = gpd.GeoDataFrame(df_lines, geometry = "geometry", crs = "EPSG:4326")

    filepath_lines = snakemake.output["lines"]
    # save substations output
    logger.info(f"Exporting clean lines to {filepath_lines}")
    parentfolder_lines = os.path.dirname(filepath_lines)
    if not os.path.exists(parentfolder_lines):
        # Create the folder and its parent directories if they don't exist
        os.makedirs(parentfolder_lines)

    gdf_lines.to_file(filepath_lines, driver="GeoJSON")
    

    ########
    ########
    ########


    fig = Figure(width = "50%", height = 600)

    m = gdf_substations.explore(name = "Buses", color = "red")
    m = gdf_lines.explore(m = m, name = "Lines")

    folium.LayerControl(collapsed = False).add_to(m)

    fig.add_child(m)
    m

    gdf_substations.explore()


    output = str(snakemake.output)
    clean_osm_data(output)