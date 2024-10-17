# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import logging
import os
import string

import itertools
import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from pyproj import Transformer
from shapely import prepare
from shapely.algorithms.polylabel import polylabel
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import linemerge, nearest_points, split
from tqdm import tqdm

logger = logging.getLogger(__name__)


GEO_CRS = "EPSG:4326"
DISTANCE_CRS = "EPSG:3035"
BUS_TOL = 500 # meters
LINES_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "circuits",
    "length",
    "underground",
    "under_construction",
    "geometry",
]
LINKS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "length",
    "underground",
    "under_construction",
    "geometry",
]
TRANSFORMERS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage_bus0",
    "voltage_bus1",
    "geometry",
]


# TODO: Add transformer s_nom capacity based on largest sum of connected line capacities
def _add_transformers(buses):
    buses_all = buses.copy().set_index("bus_id")
    logger.info("Adding unique transformers between buses of different voltage levels at each station.")

    all_transformers = gpd.GeoDataFrame(columns=TRANSFORMERS_COLUMNS+["station_id", "transformer_id"], crs=GEO_CRS)

    for g_name, g_value in buses_all.groupby("station_id"):
        if g_value["voltage"].nunique() > 1:
            combinations = list(itertools.combinations(sorted(g_value.index), 2))

            station_transformers = pd.DataFrame(
                combinations, columns=["bus0", "bus1"]
            )
            station_transformers["voltage_bus0"] = station_transformers["bus0"].map(
                g_value["voltage"]
            )
            station_transformers["voltage_bus1"] = station_transformers["bus1"].map(
                g_value["voltage"]
            )
            station_transformers["geometry"] = station_transformers.apply(
                lambda row: LineString(
                    [g_value.loc[row["bus0"]]["geometry"], g_value.loc[row["bus1"]]["geometry"]]
                ),
                axis=1,
            )
            station_transformers["station_id"] = g_name
            station_transformers["transformer_id"] = (
                g_name + "-" +
                station_transformers["voltage_bus0"].div(1e3).astype(int).astype(str) + "-" + station_transformers["voltage_bus1"].div(1e3).astype(int).astype(str)
            )
            station_transformers = gpd.GeoDataFrame(station_transformers, crs=GEO_CRS)          
            
            all_transformers = pd.concat([all_transformers, station_transformers])

    all_transformers.groupby("station_id").size()
    # Show all_transformers with more than two transformers per station
    all_transformers[all_transformers.transformer_id.isin(all_transformers.groupby("station_id").filter(lambda x: len(x) ==3)["transformer_id"])].explore()
    logger.info(f"Added {len(all_transformers)} transformers to {len(all_transformers['station_id'].unique())} stations.")

    return all_transformers[["transformer_id", "station_id"]+TRANSFORMERS_COLUMNS]


def _merge_buses_to_stations(buses, stations, distance_crs=DISTANCE_CRS, geo_crs=GEO_CRS):
    """
    
    """
    # Merge buses with same voltage and within tolerance
    logger.info(f"Merging buses of the same substation.")
    # bus types (AC != DC)
    buses_all = buses.copy().reset_index(drop=True)
    stations_all = stations.copy()
    stations_all["polygon"] = stations_all["geometry"].copy()

    # Set station ids based on station seeds (polygons)
    logger.info(" - Merging buses with the same substation id")
    buses_all = gpd.sjoin(buses_all, stations_all, how="left", predicate="within")

    # TODO: For now also include DC buses in the merging process. In the future, try to keep original HVDC converter location within substation
    buses_all = buses_all.drop_duplicates(subset=["station_id", "voltage", "country"])


    # Offsetting geometries within same substations for each voltage level
    offset = 15 # meters (radius)
    geo_to_dist = Transformer.from_crs(geo_crs, distance_crs, always_xy=True)
    dist_to_geo = Transformer.from_crs(distance_crs, geo_crs, always_xy=True)

    for g_name, g_value in tqdm(buses_all.groupby("station_id"), desc="Updating geometries", unit=" substations"):
        voltages = sorted(g_value["voltage"].unique(), reverse=True) # Sort voltags in descending order

        if len(voltages) > 1:
            poi_x, poi_y = geo_to_dist.transform(g_value["poi"].values[0].x, g_value["poi"].values[0].y)

            for idx, v in enumerate(voltages):
                poi_x_offset = poi_x + offset * np.sin(np.pi/4 + 2*np.pi*idx/len(voltages)).round(4)
                poi_y_offset = poi_y + offset * np.cos(np.pi/4 + 2*np.pi*idx/len(voltages)).round(4)

                poi_offset = Point(
                    dist_to_geo.transform(poi_x_offset, poi_y_offset)
                )

                # Update bus_name
                g_value.loc[g_value["voltage"]==v, "bus_id"] = g_name + "-" + str(int(v/1000))

                # Update geometry
                g_value.loc[g_value["voltage"]==v, "geometry"] = poi_offset

            buses_all.loc[g_value.index, "bus_id"] = g_value["bus_id"]
            buses_all.loc[g_value.index, "geometry"] = g_value["geometry"]
        else:
            v = voltages[0]
            buses_all.loc[g_value.index, "bus_id"] = g_name + "-" + str(int(v/1000))
            buses_all.loc[g_value.index, "geometry"] = g_value["poi"]

    return buses_all


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


# TODO: Last old function to improve, either vectorise or parallelise
def fix_overpassing_lines(lines, buses, distance_crs=DISTANCE_CRS, tol=1):
    """
    Fix overpassing lines by splitting them at nodes within a given tolerance,
    to include the buses being overpassed.

    Parameters:
    - lines (GeoDataFrame): The lines to be fixed.
    - buses (GeoDataFrame): The buses representing nodes.
    - distance_crs (str): The coordinate reference system (CRS) for distance calculations.
    - tol (float): The tolerance distance in meters for determining if a bus is within a line.
    Returns:
    - lines (GeoDataFrame): The fixed lines.
    - buses (GeoDataFrame): The buses representing nodes.
    """
    lines = lines.copy()
    logger.info(f"Splitting lines over overpassing nodes (Tolerance {tol} m).")
    lines_to_add = []  # list of lines to be added
    lines_to_split = []  # list of lines that have been split

    lines_epsgmod = lines.to_crs(distance_crs)
    buses_epsgmod = buses.to_crs(distance_crs)

    # set tqdm options for substation ids
    tqdm_kwargs_substation_ids = dict(
        ascii=False,
        unit=" lines",
        total=lines.shape[0],
        desc="Splitting lines",
    )

    for l in tqdm(lines.index, **tqdm_kwargs_substation_ids):
        # bus indices being within tolerance from the line
        bus_in_tol_epsg = buses_epsgmod[
            buses_epsgmod.geometry.distance(lines_epsgmod.geometry.loc[l]) <= tol
        ]

        # Get boundary points
        endpoint0 = lines_epsgmod.geometry.loc[l].boundary.geoms[0]
        endpoint1 = lines_epsgmod.geometry.loc[l].boundary.geoms[1]

        # Calculate distances
        dist_to_ep0 = bus_in_tol_epsg.geometry.distance(endpoint0)
        dist_to_ep1 = bus_in_tol_epsg.geometry.distance(endpoint1)

        # exclude endings of the lines
        bus_in_tol_epsg = bus_in_tol_epsg[(dist_to_ep0 > tol) | (dist_to_ep1 > tol)]

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
            # update name of the line if there are multiple line segments
            df_append["line_id"] = [
                str(df_append["line_id"].iloc[0])
                + (f"-{letter}" if n_geoms > 1 else "")
                for letter in string.ascii_lowercase[:n_geoms]
            ]

            lines_to_add.append(df_append)

    if not lines_to_add:
        return lines

    df_to_add = gpd.GeoDataFrame(pd.concat(lines_to_add, ignore_index=True))
    df_to_add.set_crs(lines.crs, inplace=True)
    df_to_add.set_index(lines.index[-1] + df_to_add.index, inplace=True)

    # update length
    df_to_add["length"] = df_to_add.to_crs(distance_crs).geometry.length

    # remove original lines
    lines.drop(lines_to_split, inplace=True)
    lines = df_to_add if lines.empty else pd.concat([lines, df_to_add])
    lines = gpd.GeoDataFrame(lines.reset_index(drop=True), crs=lines.crs)

    return lines


def _create_station_seeds(buses, buses_polygon, country_shapes, tol = BUS_TOL, distance_crs = DISTANCE_CRS, geo_crs=GEO_CRS):
    # Drop all buses that have bus_id starting with "way/" or "relation/" prefix
    # They are present in polygon_buffer anyway
    logger.info("Creating aggregated stations based on substation polygons")
    columns = ["bus_id", "geometry"]
    filtered_buses = buses[~buses["bus_id"].str.startswith("way/") & ~buses["bus_id"].str.startswith("relation/")]

    buses_buffer = gpd.GeoDataFrame(data = filtered_buses[columns], geometry="geometry").set_index("bus_id")
    buses_buffer["geometry"] = buses_buffer["geometry"].to_crs(distance_crs).buffer(tol).to_crs(GEO_CRS)
    buses_buffer["area"] = buses_buffer.to_crs(distance_crs).area

    buses_polygon_buffer = gpd.GeoDataFrame(data = buses_polygon[columns], geometry="geometry").set_index("bus_id")
    buses_polygon_buffer["geometry"] = buses_polygon_buffer["geometry"].to_crs(distance_crs).buffer(tol).to_crs(GEO_CRS)
    buses_polygon_buffer["area"] = buses_polygon_buffer.to_crs(distance_crs).area

    # Find the pole of inaccessibility (PoI) for each polygon, the interior point most distant from a polygon's boundary
    # Garcia-Castellanos & Lombardo 2007
    # doi.org/10.1080/14702540801897809
    buses_polygon_buffer["poi"] = buses_polygon.set_index("bus_id")["geometry"].to_crs(DISTANCE_CRS) \
        .apply(lambda polygon: polylabel(polygon, tolerance = BUS_TOL/2)) \
        .to_crs(GEO_CRS)

    buses_all_buffer = pd.concat([buses_buffer, buses_polygon_buffer])

    buses_all_agg = gpd.GeoDataFrame(geometry=[poly for poly in buses_all_buffer.union_all().geoms], crs=GEO_CRS)

    # full spatial join
    buses_all_agg = gpd.sjoin(buses_all_agg, buses_all_buffer, how="left", predicate="intersects").reset_index()
    max_area_idx = buses_all_agg.groupby('index')['area'].idxmax()
    buses_all_agg = buses_all_agg.loc[max_area_idx]
    buses_all_agg = buses_all_agg.drop(columns=["index"])
    buses_all_agg.set_index("bus_id", inplace=True)

    # Find the PoI for all rows where buses_all_agg.poi isnaway/281014575	way/281014575	POINT (19.63755 42.01247)	AL	POLYGON ((19.64054 42.01218, 19.64049 42.01196...	POLYGON ((19.64365 42.
    poi_missing = buses_all_agg["poi"].isna()
    buses_all_agg.loc[poi_missing, "poi"] = buses_all_agg.loc[poi_missing, "geometry"].to_crs(DISTANCE_CRS) \
        .apply(lambda polygon: polylabel(polygon, tolerance = BUS_TOL/2)) \
        .to_crs(GEO_CRS)

    # Update countries based on the PoI location:
    updated_country_mapping = gpd.sjoin_nearest(
        gpd.GeoDataFrame(buses_all_agg["poi"], geometry="poi", crs=geo_crs).to_crs(DISTANCE_CRS),
        country_shapes.to_crs(DISTANCE_CRS),
        how="left"
    )["name"]
    updated_country_mapping.index.name="country"

    buses_all_agg["country"] = updated_country_mapping

    # Rename rows virtual buses that are not actual substations
    buses_to_rename = buses_all_agg.loc[
        (~buses_all_agg.index.str.startswith("way/") & ~buses_all_agg.index.str.startswith("relation/")),
        ["country", "poi"]
    ]

    buses_to_rename['lat'] = buses_to_rename['poi'].apply(lambda p: p.y)
    buses_to_rename['lon'] = buses_to_rename['poi'].apply(lambda p: p.x)

    # Now sort by country, latitude (north to south), and longitude (west to east)
    buses_to_rename = buses_to_rename.sort_values(by=['country', 'lat', 'lon'], ascending=[True, False, True])
    buses_to_rename['bus_id'] = buses_to_rename.groupby('country').cumcount() + 1
    buses_to_rename['bus_id'] = buses_to_rename['country'] + buses_to_rename['bus_id'].astype(str)

    # Dict to rename virtual buses
    dict_rename = buses_to_rename["bus_id"].to_dict()

    # Rename virtual buses in buses_all_agg with dict
    buses_all_agg.rename(index=dict_rename, inplace=True)

    # extract substring before - from index
    buses_all_agg["osm_identifier"] = buses_all_agg.index.str.split("-").str[0]
    buses_all_agg.reset_index(inplace=True)
    # count how often each value in osm_identifier occurs in column
    buses_all_agg["id_occurence"] = buses_all_agg.groupby("osm_identifier")["osm_identifier"].transform("count")

    # For each row, if id_occurence =1 set bus_id = osm_identifier
    buses_all_agg.loc[buses_all_agg["id_occurence"] == 1, "bus_id"] = buses_all_agg.loc[buses_all_agg["id_occurence"] == 1, "osm_identifier"]
    buses_all_agg.set_index("bus_id", inplace=True)
    # Rename index name to station_id
    buses_all_agg.index.name = "station_id"
    buses_all_agg.drop(columns=["area", "osm_identifier", "id_occurence"], inplace=True)

    buses_all_agg["poi_perimeter"] = buses_all_agg["poi"].to_crs(distance_crs).buffer(tol/2).to_crs(geo_crs)

    return buses_all_agg


def _create_merge_mapping(lines, buses, buses_polygon):
    logger.info("Creating mapping for merging lines with same electric parameters over virtual buses.")
    lines=lines.copy()
    buses_virtual=buses.copy()

    buses_virtual = buses_virtual[buses_virtual["bus_id"].str.startswith("virtual")]
    # Drop buses_virtual that are inside buses_virtual_polygon
    bool_intersects_buses_virtual_polygon = buses_virtual.intersects(buses_polygon.union_all())
    buses_virtual = buses_virtual[~bool_intersects_buses_virtual_polygon]

    # sjoin with lines_filtered in column "connected_lines", that intersect with buses_virtual
    buses_virtual = gpd.sjoin(
        buses_virtual,
        lines[["line_id", "geometry", "voltage", "circuits"]],
        how="left",
        predicate="touches"
    )
    # Filtering, only keep where voltage_left == voltage_right
    buses_virtual = buses_virtual[buses_virtual["voltage_left"] == buses_virtual["voltage_right"]]
    # Drop voltage_right and rename voltage_left to voltage
    buses_virtual = buses_virtual.drop(columns=["voltage_right"])
    buses_virtual = buses_virtual.rename(columns={"voltage_left": "voltage"})

    # Group by bus_id, voltage, circuits and count the number of connected lines
    buses_to_remove = buses_virtual.groupby(["bus_id", "voltage"]).size().reset_index(name="count")
    buses_to_remove = buses_to_remove[buses_to_remove["count"] == 2]

    buses_to_remove = buses_virtual[buses_virtual["bus_id"].isin(buses_to_remove["bus_id"])]
    buses_to_remove = buses_to_remove.groupby(["bus_id", "circuits", "voltage"]).size().reset_index(name="count")
    # Keep only where count == 2
    buses_to_remove = buses_to_remove[buses_to_remove["count"] == 2]
    buses_to_remove = buses_virtual[buses_virtual["bus_id"].isin(buses_to_remove["bus_id"])]

    # Group by bus_id, voltage, circuits and count the number of connected lines, add column "lines" which contains the list of lines
    buses_to_remove = buses_to_remove.groupby(["bus_id", "voltage", "circuits"]) \
        .agg({
            "line_id": list,
        }).reset_index()

    unique_lines = pd.Series(itertools.chain(*buses_to_remove["line_id"])).unique()
    lines_to_merge = lines.loc[lines["line_id"].isin(unique_lines), ["line_id", "voltage", "circuits", "length", "geometry", "tag_type", "underground"]]
    lines_to_merge_dict = [(node, row.to_dict()) for node, row in lines_to_merge.set_index("line_id").iterrows()]

    # Create networkx
    G = nx.Graph()

    # Add nodes from
    # lines_to_merge
    G.add_nodes_from(lines_to_merge_dict)

    edges = [
        (
            row["line_id"][0],  # from node
            row["line_id"][1],  # to node
            {
                "bus_id": row["bus_id"],      # shared bus
            }
        )
        for _, row in buses_to_remove.iterrows()
    ]

    # Add all edges at once
    G.add_edges_from(edges)

    connected_components = nx.connected_components(G)

    # Create empty
    merged_lines = pd.DataFrame(columns=["line_id", "circuits", "voltage", "geometry", "tag_type", "underground", "contains_lines", "contains_buses"])

    # Iterate over each connected component
    for i, component in enumerate(connected_components, 1):
        # Create a subgraph for the current component
        subgraph = G.subgraph(component)

        # Extract the first node in the component
        first_node = next(iter(component))  # Get the first node (can be arbitrary)

        # Extract the attributes for the first node
        circuits = G.nodes[first_node].get('circuits', None)  # Default to None if not present
        voltage = G.nodes[first_node].get('voltage', None)    # Default to None if not present

        # Extract the geometries for all nodes in the subgraph
        geometry = [G.nodes[node].get('geometry', None) for node in subgraph.nodes()]
        geometry = linemerge(geometry)  # Merge the geometries

        # # list of all countries
        # country = ';'.join([G.nodes[node].get('country', '') for node in subgraph.nodes()])
        # country = ';'.join(sorted(set(country.split(';'))))

        # Contains lines
        contains_lines = list(subgraph.nodes())
        number_of_lines = int(len(contains_lines))
        # Find the node with the highest "length" parameter
        node_longest = max(subgraph.nodes(), key=lambda node: G.nodes[node].get('length', 0))
        tag_type = G.nodes[node_longest].get('tag_type', None)
        underground = G.nodes[node_longest].get('underground', None)

        # Extract the list of edges (lines) in the subgraph
        contains_buses = list()
        for edge in subgraph.edges():
            contains_buses.append(G.edges[edge].get('bus_id', None))

        subgraph_data = {
            "line_id": ["merged_"+str(node_longest)+"+"+str(number_of_lines-1)],  # line_id
            "circuits": [circuits],
            "voltage": [voltage],
            "geometry": [geometry],
            "tag_type": [tag_type],
            "underground": [underground],
            "contains_lines": [contains_lines],
            "contains_buses": [contains_buses],
        }

         # Convert to DataFrame and append to the merged_lines DataFrame
        merged_lines = pd.concat([merged_lines, pd.DataFrame(subgraph_data)], ignore_index=True)
        merged_lines = gpd.GeoDataFrame(merged_lines, geometry="geometry", crs=GEO_CRS)

        # Drop all closed linestrings (circles)
        merged_lines = merged_lines[merged_lines["geometry"].apply(lambda x: not x.is_closed)]

    return merged_lines


def _merge_lines_over_virtual_buses(lines, buses, merged_lines_map):
    lines_merged = lines.copy()
    buses_merged = buses.copy()

    lines_to_remove = merged_lines_map["contains_lines"].explode().unique()
    buses_to_remove = merged_lines_map["contains_buses"].explode().unique()

    logger.info(f"Merging {len(lines_to_remove)} lines over virtual buses and dropping {len(buses_to_remove)} virtual buses.")

    # Remove lines and buses
    lines_merged = lines_merged[~lines_merged["line_id"].isin(lines_to_remove)]
    lines_merged["contains_lines"] = lines_merged["line_id"].apply(lambda x: [x])
    buses_merged = buses_merged[~buses_merged["bus_id"].isin(buses_to_remove)]
    # Add new lines from merged_lines_map to lines_merged

    ### Update lines to add
    lines_to_add = merged_lines_map.copy().reset_index(drop=True)

    # Update columns
    lines_to_add["under_construction"] = False
    lines_to_add["tag_frequency"] = 50
    lines_to_add["dc"] = False
    lines_to_add["bus0"] = None
    lines_to_add["bus1"] = None
    lines_to_add["length"] = lines_to_add["geometry"].to_crs(DISTANCE_CRS).length

    # Reorder
    lines_to_add = lines_to_add[lines_merged.columns]

    # Concatenate
    lines_merged = pd.concat([lines_merged, lines_to_add], ignore_index=True)

    return lines_merged, buses_merged


def _identify_linestring_between_polygons(multiline, polygon0, polygon1, geo_crs=GEO_CRS, distance_crs=DISTANCE_CRS):
    list_lines = [line for line in multiline.geoms] if multiline.geom_type == 'MultiLineString' else [multiline]
    prepare(polygon0)
    prepare(polygon1)
    gdf_polygon0 = gpd.GeoDataFrame(geometry=[polygon0], crs=geo_crs).to_crs(distance_crs)
    gdf_polygon1 = gpd.GeoDataFrame(geometry=[polygon1], crs=geo_crs).to_crs(distance_crs)

    idx = None
    for line in list_lines:
        prepare(line)
        gdf_line = gpd.GeoDataFrame(geometry=[line], crs=geo_crs).to_crs(distance_crs)
        touches = gdf_line.intersects(gdf_polygon0.buffer(1e-2)) & gdf_line.intersects(gdf_polygon1.buffer(1e-2))
        if touches.any():
            idx = list_lines.index(line)
            break
    
    return list_lines[idx] if idx is not None else multiline


def _remove_loops_from_multiline(multiline):
    elements_initial = [line for line in multiline.geoms] if multiline.geom_type == 'MultiLineString' else [multiline]

    if not any([line.is_closed for line in elements_initial]):
        return multiline
    
    elements = elements_initial
    iteration_count = 0
    while any([line.is_closed for line in elements]) and iteration_count < 5:
        elements = [line for line in elements if not line.is_closed]
        geometry_updated = linemerge(elements)

        elements = [line for line in geometry_updated.geoms] if geometry_updated.geom_type == 'MultiLineString' else [geometry_updated]
        iteration_count += 1

        if not any([line.is_closed for line in elements]):
            break
    
    return geometry_updated


def _map_endpoints_to_buses(connection, buses, shape="station_polygon", id_col="line_id", sjoin="intersects", distance_crs=DISTANCE_CRS, geo_crs=GEO_CRS):
    buses_all = buses.copy().set_index("bus_id")
    buses_all["station_polygon"] = buses_all["polygon"].copy()
    buses_all = gpd.GeoDataFrame(buses_all, geometry="polygon", crs = buses.crs)

    lines_all = connection.copy().set_index(id_col)

    for coord in range(2):
        # Obtain endpoints
        endpoints = lines_all[["voltage", "geometry"]].copy()
        endpoints["geometry"] = endpoints["geometry"].apply(lambda x: x.boundary.geoms[coord])
        if sjoin=="intersects":
            endpoints = gpd.sjoin(endpoints, buses_all, how="left", predicate="intersects")
            endpoints = endpoints[endpoints["voltage_left"] == endpoints["voltage_right"]]
        if sjoin=="nearest":
            endpoints = gpd.sjoin_nearest(endpoints.to_crs(distance_crs), buses_all.to_crs(distance_crs), how="left")
            # Groupby index and keep the row with highest voltage_right
            endpoints.reset_index(inplace=True)
            endpoints = endpoints.loc[endpoints.groupby("link_id")["voltage_right"].idxmax()]
            endpoints.set_index("link_id", inplace=True)

        # rename voltage_left
        endpoints = endpoints.drop(columns=["voltage_right"])
        endpoints = endpoints.rename(columns={"voltage_left": "voltage"})

        lines_all[f"poi_perimeter{coord}"] = endpoints["poi_perimeter"]
        lines_all[f"station_polygon{coord}"] = endpoints["station_polygon"]
        lines_all[f"bus{coord}"] = endpoints["bus_id"]
    
    lines_all["geometry"] = lines_all.apply(lambda row: row["geometry"].difference(row[shape+"0"]), axis=1)
    lines_all["geometry"] = lines_all.apply(lambda row: row["geometry"].difference(row[shape+"1"]), axis=1)

    # Drop lines that have same bus0 and bus1
    lines_all = lines_all[lines_all["bus0"] != lines_all["bus1"]]

    contains_stubs = lines_all["geometry"].apply(lambda x: isinstance(x, MultiLineString))
    if contains_stubs.any():
        lines_stubs = lines_all.loc[contains_stubs].copy()
        lines_stubs["linestrings"] = lines_stubs["geometry"].apply(lambda x: [line for line in x.geoms] if x.geom_type == 'MultiLineString' else [x])

        # Remove closed subgeometries (circles), if any
        lines_stubs["geometry"] = lines_stubs["geometry"].apply(_remove_loops_from_multiline)

        if shape == "station_polygon":
            # Only keep linestrings that are actually connection shape of bus0 and bus1
            lines_stubs["geometry"] = lines_stubs.apply(lambda row:
                _identify_linestring_between_polygons(
                    row["geometry"], 
                    row[f"{shape}0"], 
                    row[f"{shape}1"],
                ),
                axis=1
            )
            lines_all.loc[lines_stubs.index, "geometry"] = lines_stubs["geometry"]
                                                    
        # If the cutton should be done at the poi_perimeters
        if shape == "poi_perimeter":
            for coord in range(2):
                lines_stubs = lines_all.loc[contains_stubs].copy()
                lines_stubs["linestrings"] = lines_stubs["geometry"].apply(lambda x: [line for line in x.geoms] if x.geom_type == 'MultiLineString' else [x])

                # Check for each individual element in list of linestrings, if they are within the station_polygon, if yes delete
                lines_stubs["linestrings"] = lines_stubs.apply(lambda row: [line for line in row["linestrings"] if not row[f"station_polygon{coord}"].contains(line)], axis=1)

                # Update geometry through linemerge
                lines_stubs["geometry"] = lines_stubs["linestrings"].apply(lambda x: linemerge(x))

                # Update lines_all
                lines_all.loc[lines_stubs.index, "geometry"] = lines_stubs["geometry"]

    lines_all.reset_index(inplace=True)
    lines_all.drop(columns=["poi_perimeter0", "poi_perimeter1", "station_polygon0", "station_polygon1"], inplace=True)

    return lines_all


def _add_point_to_line(linestring, point):
    """
    Adds the bus coordinate to a linestring by extending the linestring with a
    new segment.

    Parameters:
    linestring (LineString): The original linestring to extend.
    point (Point): The shapely.Point of the bus.

    Returns:
    merged (LineString): The extended linestring with the new segment.
    """
    start = linestring.boundary.geoms[0]
    end = linestring.boundary.geoms[1]

    dist_to_start = point.distance(start)
    dist_to_end = point.distance(end)

    if dist_to_start < dist_to_end:
        new_segment = LineString([point, start])
    else:
        new_segment = LineString([point, end])

    merged = linemerge([linestring, new_segment])

    return merged


def _extend_lines_to_buses(connection, buses):
    lines_all = connection.copy()
    buses_all = buses.copy()

    logger.info("Extending line/link geometry to mapped bus points.")
    lines_all = lines_all.merge(buses_all[["geometry", "bus_id"]], left_on="bus0", right_on="bus_id", how="left", suffixes=("", "_right"))
    lines_all.drop(columns=["bus_id"], inplace=True)
    lines_all.rename(columns={"geometry_right": "bus0_point"}, inplace=True)

    lines_all = lines_all.merge(buses_all[["geometry", "bus_id"]], left_on="bus1", right_on="bus_id", how="left", suffixes=("", "_right"))
    lines_all.drop(columns=["bus_id"], inplace=True)
    lines_all.rename(columns={"geometry_right": "bus1_point"}, inplace=True)

    lines_all["geometry"] = lines_all.apply(lambda row: _add_point_to_line(row["geometry"], row["bus0_point"]), axis=1)
    lines_all["geometry"] = lines_all.apply(lambda row: _add_point_to_line(row["geometry"], row["bus1_point"]), axis=1)

    # Drop bus0_point and bus1_point
    lines_all.drop(columns=["bus0_point", "bus1_point"], inplace=True)

    return lines_all


def _merge_identical_lines(lines):
    lines_all = lines.copy()
    lines_to_drop = []

    logger.info("Aggregating lines with identical geometries and voltage levels.")
    for g_name, g_value in tqdm(lines_all.groupby(["geometry", "voltage"]), desc="Aggregating identical lines", unit=" lines"):
        line_ids = list(g_value["line_id"])

        if len(line_ids) > 1:
            # New aggregated line
            lid_old = line_ids[0]
            lid_agg = lid_old.split("-")[0]
            circuits_agg = g_value["circuits"].sum()

            # Update line with aggregated parameters
            lines_all.loc[lines_all["line_id"] == lid_old, "line_id"] = lid_agg
            lines_all.loc[lines_all["line_id"] == lid_agg, "circuits"] = circuits_agg

            # Add other lines to lines to drop
            lines_to_drop += line_ids[1:]

    logger.info(f"In total {len(lines_to_drop)} lines aggregated.")

    # Drop lines
    lines_all = lines_all[~lines_all["line_id"].isin(lines_to_drop)]

    # Update line ids to make them unique again
    # Add voltage suffix to line_id
    lines_all["line_id"] = lines_all["line_id"] + "-" + lines_all["voltage"].div(1e3).astype(int).astype(str)

    return lines_all


def _add_line_endings(buses, lines, add=0, name="line-end"):
    buses_all = buses.copy()

    endpoints0 = lines[["voltage", "geometry"]].copy()
    endpoints0["geometry"] = endpoints0["geometry"].apply(lambda x: x.boundary.geoms[0])

    endpoints1 = lines[["voltage", "geometry"]].copy()
    endpoints1["geometry"] = endpoints1["geometry"].apply(lambda x: x.boundary.geoms[1])

    endpoints = pd.concat([endpoints0, endpoints1], ignore_index=True)
    endpoints.drop_duplicates(subset=["geometry", "voltage"], inplace=True)
    endpoints.reset_index(drop=True, inplace=True)


    endpoints["bus_id"] = endpoints.index+add+1
    endpoints["bus_id"] = "virtual" + "-" + endpoints["bus_id"].astype(str)

    endpoints["tag_source"] = name

    return endpoints[["bus_id", "voltage", "geometry", "tag_source"]]


def _map_links_to_dc_buses(links, dc_buses):
    links_all = links.copy().set_index("link_id")
    max_distance = 120000 # m, arbitrary maximum relevant for islands, and if second country is missing
    dc_buses_all=dc_buses.copy().set_index("bus_id")
    dc_buses_polygon = gpd.GeoDataFrame(dc_buses_all["polygon"], geometry="polygon", crs=dc_buses.crs)

    links_all["bus0"] = links_all["geometry"].apply(lambda x: x.boundary.geoms[0])
    links_all["bus1"] = links_all["geometry"].apply(lambda x: x.boundary.geoms[1])

    links_all["bus0"] = gpd.sjoin_nearest(
        gpd.GeoDataFrame(links_all["bus0"], geometry="bus0", crs=links_all.crs).to_crs(DISTANCE_CRS),
        dc_buses_polygon.to_crs(DISTANCE_CRS),
        how="left",
        max_distance=max_distance, 
    )["bus_id"]

    links_all["bus1"] = gpd.sjoin_nearest(
        gpd.GeoDataFrame(links_all["bus1"], geometry="bus1", crs=links_all.crs).to_crs(DISTANCE_CRS),
        dc_buses_polygon.to_crs(DISTANCE_CRS),
        how="left",
        max_distance=max_distance, 
    )["bus_id"]

    # list of dc_buses that have been used
    mapped_dc_buses = links_all[["bus0", "bus1"]].stack().unique()
    # Drop remaining dc_buses
    dc_buses_all = dc_buses_all[dc_buses_all.index.isin(mapped_dc_buses)]

    # Obtain voltage for DC buses
    bus_voltages0 = links_all.set_index("bus0")["voltage"]
    bus_voltages1 = links_all.set_index("bus1")["voltage"]
    bus_voltages = pd.DataFrame(pd.concat([bus_voltages0, bus_voltages1], axis=0))
    bus_voltages = bus_voltages.groupby(bus_voltages.index).max()
    dc_buses_all["voltage"] = bus_voltages

    # Reset index
    links_all.reset_index(inplace=True)
    dc_buses_all.reset_index(inplace=True)

    return links_all, dc_buses_all


def _add_dc_buses(converters_polygon, links, buses, distance_crs=DISTANCE_CRS, geo_crs=GEO_CRS):
    dc_buses = converters_polygon.copy()
    max_distance = 5000

    # Rough filtering: Only keep dc_buses if distance to links union is less than 5000 meters
    links_closeby = dc_buses.to_crs(distance_crs).distance(links.to_crs(distance_crs).union_all()) < 5000
    dc_buses = dc_buses[links_closeby]

    dc_buses.rename(columns={"id": "bus_id", "geometry": "polygon"}, inplace=True)

    # Determine PoI for each dc bus
    dc_buses["poi"] = dc_buses["polygon"].to_crs(distance_crs) \
        .apply(lambda polygon: polylabel(polygon, tolerance = BUS_TOL/2)) \
        .to_crs(geo_crs)
    dc_buses["geometry"] = dc_buses["poi"]

    # Map dc_buses to stations 
    ac_buses = buses.copy().set_index("bus_id")
    # Group ac_buses by station, keep row with highest voltage
    ac_buses = ac_buses.loc[ac_buses.groupby("station_id")["voltage"].idxmax()]
    ac_buses_polygon = gpd.GeoDataFrame(ac_buses[["station_id", "polygon"]], geometry="polygon", crs=ac_buses.crs)

    dc_buses = gpd.sjoin_nearest(
        dc_buses.to_crs(distance_crs),
        ac_buses_polygon.to_crs(distance_crs),
        how="left",
        lsuffix=None,
        rsuffix="ac",
        max_distance=max_distance,
    ).to_crs(geo_crs)

    return dc_buses

def _add_converter_links(dc_buses, buses):
    converter_links = dc_buses.copy()
    # Drop rows that do not have a bus_id_ac
    converter_links = converter_links[~converter_links["bus_id_ac"].isna()]

    ac_buses_geometry = buses.copy().set_index("bus_id")[["geometry"]]
    ac_buses_geometry = ac_buses_geometry[ac_buses_geometry.index.isin(converter_links["bus_id_ac"])]

    converter_links.rename(columns={"bus_id": "converter_id"}, inplace=True)
    converter_links["bus0"] = converter_links["converter_id"]
    converter_links["bus1"] = converter_links["bus_id_ac"]
    converter_links["ac_bus_geometry"] = converter_links["bus_id_ac"].apply(lambda x: ac_buses_geometry.loc[x])

    converter_links["geometry"] = converter_links.apply(
        lambda row: LineString([row["geometry"], row["ac_bus_geometry"]]),
        axis=1
    )

    converter_links["converter_id"] = "conv-" + converter_links["converter_id"] 

    return converter_links[["converter_id", "bus0", "bus1", "voltage", "geometry"]]


def build_network(inputs, countries):
    logger.info("Reading input data.")
    buses = gpd.read_file(inputs["substations"])
    buses = buses[["bus_id", "voltage", "geometry", "tag_source"]]
    buses_polygon = gpd.read_file(inputs["substations_polygon"])

    # Simplify buses_polygon
    buses_polygon["bus_id"] = buses_polygon["bus_id"].apply(lambda x: x.split("-")[0])
    buses_polygon.drop_duplicates(subset=["bus_id", "country", "geometry"], inplace=True)
    buses_polygon.drop(columns=["voltage"], inplace=True)

    lines = gpd.read_file(inputs["lines"])
    lines = lines.drop(columns=["dc", "country"])

    lines = _merge_identical_lines(lines)

    # Add line endings
    # Drop all buses that do not intersect with lines
    buses = buses[buses.intersects(lines.union_all())]
    buses_line_endings = _add_line_endings(buses, lines)

    buses = pd.concat([buses, buses_line_endings], ignore_index=True)

    # Split lines overpassing nearby buses (tolerance 1 m)
    lines = fix_overpassing_lines(lines, buses, DISTANCE_CRS, tol=1)

    # Update end points
    bool_virtual_buses = buses["bus_id"].str.startswith("virtual")
    buses = buses[~bool_virtual_buses]
    buses_updated_line_endings = _add_line_endings(buses, lines)
    buses = pd.concat([buses, buses_updated_line_endings], ignore_index=True)

    merged_lines_map = _create_merge_mapping(lines, buses, buses_polygon)
    lines, buses = _merge_lines_over_virtual_buses(lines, buses, merged_lines_map)

    # Create station seeds
    stations = _create_station_seeds(buses, buses_polygon, country_shapes)
    buses = _merge_buses_to_stations(buses, stations)

    # Drop lines that are fully within stations.polygon
    internal_lines = gpd.sjoin(lines, stations, how="inner", predicate="within").line_id
    logger.info(f"Dropping {len(internal_lines)} lines that are fully within aggregated substations.")
    lines = lines[~lines.line_id.isin(internal_lines)].reset_index(drop=True)
    logger.info(f"Keeping {len(lines)} lines.")

    lines = _map_endpoints_to_buses(
        connection=lines,
        buses=buses,
        shape="station_polygon",
        id_col="line_id",
        sjoin="intersects"
    )
    lines = _extend_lines_to_buses(lines, buses) # Extend lines

    # Update all lengths
    lines["length"] = lines.to_crs(DISTANCE_CRS).length

    # Add transformers
    transformers = _add_transformers(buses)

    ### Building steps for HVDC links
    # Debug links
    links = gpd.read_file(inputs["links"])
    converters_polygon = gpd.read_file(inputs["converters_polygon"])

    # Create DC buses
    dc_buses = _add_dc_buses(converters_polygon, links, buses)
    links, dc_buses = _map_links_to_dc_buses(links, dc_buses)

    # Add suffix
    links["link_id"] = links["link_id"] + "-" + links["voltage"].div(1e3).astype(int).astype(str) + "-DC"

    # Extend DC links to DC buses
    links = _extend_lines_to_buses(links, dc_buses)

    # Add PyPSA converter links between DC buses and AC buses
    links_converters = _add_converter_links(dc_buses, buses)
    
    # TODO list of accepted voltages per country

    ### Debugging 
    # TODO: Delete
    map = None
    map = stations.loc[(stations.index.str.startswith("way") | stations.index.str.startswith("relation"))].explore(color="red")
    map = stations.loc[~(stations.index.str.startswith("way") | stations.index.str.startswith("relation"))].explore(color="green", m=map)
    map = buses_polygon.explore(m=map, color="yellow", popup=True)
    map = lines.explore(m=map, color="blue",popup=True)
    map = links.explore(m=map, color="purple", popup=True)
    map = buses.explore(m=map, color="red", popup=True)
    # map = stations.poi.explore(m=map, color="red", popup=True)
    map = transformers.explore(m=map, color="orange", popup=True)
    map = dc_buses.polygon.explore(m=map, color="black")
    map = dc_buses.explode().explore(m=map, color="black", popup=True)
    map = links_converters.explore(m=map, color="pink", popup=True)
    map

    logger.info("Saving outputs")

    ### Convert output to pypsa-eur friendly format
    # Rename "substation" in buses["symbol"] to "Substation"
    buses["symbol"] = buses["symbol"].replace({"substation": "Substation"})

    # Drop unnecessary index column and set respective element ids as index
    lines.set_index("line_id", inplace=True)
    if not links.empty:
        links.set_index("link_id", inplace=True)
    converters.set_index("converter_id", inplace=True)
    transformers.set_index("transformer_id", inplace=True)
    buses.set_index("bus_id", inplace=True)

    # Convert voltages from V to kV
    lines["voltage"] = lines["voltage"] / 1000
    if not links.empty:
        links["voltage"] = links["voltage"] / 1000
    if not converters.empty:
        converters["voltage"] = converters["voltage"] / 1000
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
    lines = lines[LINES_COLUMNS]
    if not links.empty:
        links = links[LINKS_COLUMNS]
    else:
        links = pd.DataFrame(columns=["link_id"] + LINKS_COLUMNS)
        links.set_index("link_id", inplace=True)
    transformers = transformers[TRANSFORMERS_COLUMNS]

    return buses, converters, lines, links, transformers


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_osm_network")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    countries = snakemake.params.countries
    country_shapes = gpd.read_file(snakemake.input["country_shapes"]).set_index("name")

    buses, converters, lines, links, transformers = build_network(
        snakemake.input,
        countries,
    )

    # Export to csv for base_network
    buses.to_csv(snakemake.output["substations"], quotechar="'")
    lines.to_csv(snakemake.output["lines"], quotechar="'")
    links.to_csv(snakemake.output["links"], quotechar="'")
    converters.to_csv(snakemake.output["converters"], quotechar="'")
    transformers.to_csv(snakemake.output["transformers"], quotechar="'")

    # Export to GeoJSON for quick validations
    buses.to_file(snakemake.output["substations_geojson"])
    lines.to_file(snakemake.output["lines_geojson"])
    links.to_file(snakemake.output["links_geojson"])
    converters.to_file(snakemake.output["converters_geojson"])
    transformers.to_file(snakemake.output["transformers_geojson"])
