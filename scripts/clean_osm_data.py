# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
TODO To fill later
"""

import geopandas as gpd
import json
import logging
import pandas as pd
import re
from shapely.geometry import LineString, Point
import tqdm.auto as tqdm

from _helpers import configure_logging
logger = logging.getLogger(__name__)

def clean_osm_data(output):
    with open(output, "w") as file:
        file.write("Hello, world!\n")


def _create_linestring(row):
    coords = [(coord['lon'], coord['lat']) for coord in row["geometry"]]
    return LineString(coords)


def _clean_voltage(column):
    """
    Function to clean the raw voltage column: manual fixing and drop nan values

    Args:
    - column: pandas Series, the column to be cleaned

    Returns:
    - column: pandas Series, the cleaned column
    """
    column = (
        column
        .astype(str)
        .str.lower()
        .str.replace("fixme", "")
        .str.replace("(temp 150000)", "")
        .str.replace("low", "1000")
        .str.replace("minor", "1000")
        .str.replace("medium", "33000")
        .str.replace("med", "33000")
        .str.replace("m", "33000")
        .str.replace("high", "150000")
        .str.replace("unknown", "")
        .str.replace("23000-109000", "109000")
        .str.replace("INF", "")
        .str.replace("<", "")
        .str.replace("?", "")
        .str.replace(",", "")
        .str.replace(" ", "")
        .str.replace("_", "")
        .str.replace("kv", "000")
        .str.replace("v", "")
        .str.replace("/", ";") 
        .str.replace("nan", "")
        .str.replace("<NA>", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r'[^0-9;]', '', x))

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
        .str.replace("1/3", "1")
        .str.replace("<NA>", "")
        .str.replace("nan", "")
    )

    # Remove all remaining non-numeric characters except for semicolons
    column = column.apply(lambda x: re.sub(r'[^0-9;]', '', x))

    column.dropna(inplace=True)
    return column.astype(str)


def _clean_frequency(column):
    column = column.copy()
    to_fifty = column.astype(str) != "0"
    column[to_fifty] = "50"    

    return column


def _split_voltage(df):
    to_split = df['voltage'].str.contains(';')
    new_rows = []
    for index, row in df[to_split].iterrows():
        split_values = row["voltage"].split(';')
        new_sub_id_len = int(len(split_values))
        for i, value in enumerate(split_values):
            new_sub_id = str(i+1)
            new_id = str(row['id']) + '_' + new_sub_id
            new_row = {
                'id': new_id, 
                'sub_id': new_sub_id,
                'sub_id_len': new_sub_id_len,
                'bounds': row['bounds'],
                'nodes': row['nodes'],
                'geometry': row['geometry'],
                'power': row['power'],
                'cables': row['cables'],
                'circuits': row['circuits'],
                'frequency': row['frequency'],
                'voltage': value, 
                'wires': row['wires'],}
            new_rows.append(new_row)

    # Create DataFrame from split rows
    split_df = pd.DataFrame(new_rows)
    df_new = pd.concat([df[~to_split], split_df])
    df_new["sub_id_len"] = df_new["sub_id_len"].astype(int)

    # Append the original DataFrame with split_df
    return df_new


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_osm_data")
    
    configure_logging(snakemake)
    logger.info("Dummy log: clean_osm_data()")

    # input_path = snakemake.input.lines_way + snakemake.input.cables_way
    # input_path = {
    #     "lines": snakemake.input.lines_way,
    #     "cables": snakemake.input.cables_way,
    # }

    # columns = ["id", "sub_id", "sub_id_len", "bounds", "nodes", "geometry", "power", "cables", "circuits", "frequency", "voltage", "wires"]
    # df_lines = pd.DataFrame(columns=columns)
    # crs = "EPSG:4326"

    # # using tqdm loop over input path

    # for key in input_path:
    #     logger.info(f"Processing {key}...")
    #     for idx, ip in enumerate(input_path[key]):
    #         if os.path.exists(ip) and os.path.getsize(ip) > 400: # unpopulated OSM json is about 51 bytes
    #             logger.info(f" - Importing {key} {str(idx+1).zfill(2)}/{str(len(input_path[key])).zfill(2)}: {ip}")
    #             with open(ip, "r") as f:
    #                 data = json.load(f)
                
    #             df = pd.DataFrame(data['elements'])
    #             df["id"] = df["id"].astype(str)
    #             df["sub_id"] = "0" # initiate sub_id column with 0
    #             df["sub_id_len"] = 0 # initiate sub_id column with 0

    #             col_tags = ["power", "cables", "circuits", "frequency", "voltage", "wires"]

    #             tags = pd.json_normalize(df["tags"]) \
    #                 .map(lambda x: str(x) if pd.notnull(x) else x)
                
    #             for ct in col_tags:
    #                 if ct not in tags.columns:
    #                     tags[ct] = pd.NA
                
    #             tags = tags.loc[:, col_tags]

    #             df = pd.concat([df, tags], axis="columns") 
    #             df.drop(columns=["type", "tags"], inplace=True)
                
    #             df_lines = pd.concat([df_lines, df], axis="rows")

    #         else:
    #             logger.info(f" - Skipping {key} {str(idx+1).zfill(2)}/{str(len(input_path[key])).zfill(2)} (empty): {ip}")
    #             continue
    #     logger.info("---")
    
    # # Drop duplicates
    # df_lines.drop_duplicates(subset="id", inplace=True)

    # df_lines["voltage"] = _clean_voltage(df_lines["voltage"])
    # # drop voltage = ""
    # df_lines = _split_voltage(df_lines)
    # df_lines = df_lines[df_lines["voltage"] != ""]
    # df_lines["voltage"] = df_lines["voltage"].astype(int, errors="ignore")

    # # Drop voltages below 220 kV
    # df_lines = df_lines[df_lines["voltage"] >= 220000]

    # # Clean frequencies
    # df_lines["frequency"] = _clean_frequency(df_lines["frequency"])
    # df_lines["frequency"] = df_lines["frequency"].astype(int, errors="ignore")

    # # Clean circuits
    # df_lines["circuits"] = _clean_circuits(df_lines["circuits"])
    # # Map correct circuits to lines that where split
    
    # # Initiate new column for cleaned circuits with values that are already valid:
    # # Condition 1: Length of sub_id is 0, the line was not split
    # # Condition 2: Number of entries in circuits separated by semicolon is 1, value is unique
    # # Condition 3: Circuits is not an empty string
    # # Condition 4: Circuits is not "0"
    # bool_circuits_valid = (df_lines["sub_id_len"] == 0) & \
    #     (df_lines["circuits"].apply(lambda x: len(x.split(";"))) == 1) & \
    #     (df_lines["circuits"] != "") & \
    #     (df_lines["circuits"] != "0")
        
    # df_lines.loc[bool_circuits_valid, "circuits_clean"] = df_lines.loc[bool_circuits_valid, "circuits"]
    
    # # Boolean to check if sub_id_len is equal to the number of circuits
    # bool_equal = df_lines["sub_id_len"] == df_lines["circuits"] \
    #                 .apply(lambda x: len(x.split(";")))
    # op_equal = lambda row: row["circuits"].split(";")[int(row["sub_id"])-1]
        
    # df_lines.loc[bool_equal, "circuits_clean"] = df_lines[bool_equal] \
    #     .apply(op_equal, axis=1)
    
    # bool_larger = df_lines["sub_id_len"] > \
    #     df_lines["circuits"].apply(lambda x: len(x.split(";")))
    
    # pd.set_option('display.max_rows', None)
    # df_lines.loc[bool_larger, ["id", "sub_id", "sub_id_len", "cables", "circuits", "circuits_clean", "frequency"]]





    # df_lines[df_lines["sub_id_len"] > 0]["circuits"]


    # df_lines["geometry"] = df_lines.apply(_create_linestring, axis=1)    
    # gdf = gpd.GeoDataFrame(
    #     df_lines[["id", "sub_id", "sub_id_len", "power", "cables", "circuits", "voltage", "geometry"]], 
    #     geometry = "geometry", crs = "EPSG:4326"
    #     )
    
    # gdf.explore()
    # df_lines.voltage.unique()

    # df_lines.circuits.apply(lambda x: x.split(";")).explode().unique()

    # ol_lines_way = ["id", "power", "cables", "circuits", "frequency", "voltage"]

    # # gdf = gpd.read_file(lines_way[3])
    # # gdf2 = gpd.GeoDataFrame(gdf, geometry=gdf.geometry)
    # # df = gdf.to_json()

    # # gdf.to_file("example.geojson", layer_options={"ID_GENERATE": "YES"})


    output = str(snakemake.output)
    clean_osm_data(output)




# # Example DataFrame
# data = {'id': ["ID1", "ID2", "ID3", "ID4", "ID5"],
#         'A': ["220000", "380000", ";100000", "220000;220000;380000", "220000;;400000;700000"],
#         'B': [1, 2, 3, 4, 5],
#         'C': [6, 7, 8, 9, 10]}
# df = pd.DataFrame(data)

# # Split the entries in column A that contain a semicolon
# split_rows = df[df['A'].str.contains(';')]
# split_values = split_rows['A'].str.split(';', expand=True)

# # Create two copies of the rows containing semicolons, one for each split value
# split_rows_1 = split_rows.copy()
# split_rows_2 = split_rows.copy()

# # Update column A in the split rows to contain the split values
# split_rows_1['A'] = split_values[0]
# split_rows_2['A'] = split_values[1]

# # Concatenate the split rows with the original DataFrame, excluding the rows containing semicolons
# result_df = pd.concat([df[~df.index.isin(split_rows.index)], split_rows_1, split_rows_2], ignore_index=True)

# # Display the result
# print(result_df)


# '# Sample DataFrame
# data = {'id': ["ID1", "ID2", "ID3", "ID4", "ID5"],
#         'voltage': ["220000", "380000", ";100000", "220000;220000;380000", "220000;;400000;700000"],
#         'B': [1, 2, 3, 4, 5],
#         'C': [6, 7, 8, 9, 10]}
# df = pd.DataFrame(data)

# # Find rows to split
# to_split = df['voltage'].str.contains(';')

# # Splitting entries and creating new rows


# new_rows = []

# for index, row in df[to_split].iterrows():
#     split_values = row["voltage"].split(';')
#     for i, value in enumerate(split_values):
#         new_id = str(row['id']) + '_' + str(i+1)
#         new_row = {
#             'id': new_id, 
#             'bounds': row['bounds'],
#             'nodes': row['nodes'],
#             'geometry': row['geometry'],
#             'cables': row['cables'],
#             'circuits': row['circuits'],
#             'frequency': row['frequency'],
#             'voltage': value, 
#             'wires': row['wires'],}
#         new_rows.append(new_row)

# # Create DataFrame from split rows
# split_df = pd.DataFrame(new_rows)

# # Append the original DataFrame with split_df
# final_df = pd.concat([df[~to_split], split_df])

# print(final_df)



# from shapely.geometry import LineString
# import numpy as np
# import matplotlib.pyplot as plt

# def offset_line(original_line, distance):
#     # Compute the direction vector between the two endpoints
#     direction_vector = np.array(original_line.coords[1]) - np.array(original_line.coords[0])

#     # Compute the orthogonal vector
#     orthogonal_vector = np.array([-direction_vector[1], direction_vector[0]])

#     # Normalize the orthogonal vector
#     orthogonal_vector /= np.linalg.norm(orthogonal_vector)

#     # Compute the offset LineString
#     offset_points = []
#     for point in original_line.coords:
#         offset_point = np.array(point) + distance * orthogonal_vector
#         offset_points.append((offset_point[0], offset_point[1]))

#     return LineString(offset_points)

# # Example usage:
# original_line = lines.iloc[5]
# offset_distance = 1.0
# b = offset_line(original_line, offset_distance)

# # Plot both LineStrings
# fig, ax = plt.subplots()
# x, y = original_line.xy
# ax.plot(x, y, label='Original LineString')
# x, y = offset_line.xy
# ax.plot(x, y, label='Offset LineString')
# ax.set_aspect('equal')
# ax.legend()
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('Original and Offset LineStrings')
# plt.grid(True)
# plt.show()