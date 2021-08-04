"""
Builds clustered natural gas network based on data from:

    [1] the SciGRID Gas project
        (https://www.gas.scigrid.de/)

    [2] ENTSOG capacity map
        (https://www.entsog.eu/sites/default/files/2019-10/Capacities%20for%20Transmission%20Capacity%20Map%20RTS008_NS%20-%20DWH_final.xlsx)
"""

import logging
logger = logging.getLogger(__name__)

import re
import json

import pandas as pd
import geopandas as gpd
import numpy as np

from shapely.geometry import Point


def concat_gdf(gdf_list, crs='EPSG:4326'):
    """Convert to gepandas dataframe with given Coordinate Reference System (crs)."""
    return gpd.GeoDataFrame(pd.concat(gdf_list),crs=crs)


def string2list(string, with_None=True):
    """Convert string format to a list."""
    p = re.compile('(?<!\\\\)\'')
    string = p.sub('\"', string)

    if with_None:
        p2 = re.compile('None')
        string = p2.sub('\"None\"', string)

    return json.loads(string)


def load_gas_network(df_path):
    """Load and format gas network data."""

    df = pd.read_csv(df_path, sep=',')

    df.long = df.long.apply(string2list)
    df.lat = df.lat.apply(string2list)
    df.node_id = df.node_id.apply(string2list)

    # pipes which can be used in both directions
    both_direct_df = df[df.is_bothDirection == 1].reset_index(drop=True)
    both_direct_df.node_id = both_direct_df.node_id.apply(lambda x: [x[1], x[0]])
    both_direct_df.long = both_direct_df.long.apply(lambda x: [x[1], x[0]])
    both_direct_df.lat = both_direct_df.lat.apply(lambda x: [x[1], x[0]])

    df_singledirect = pd.concat([df, both_direct_df]).reset_index(drop=True)
    df_singledirect.drop('is_bothDirection', axis=1)

    # create shapely geometry points
    df['point1'] = df.apply(lambda x: Point((x['long'][0], x['lat'][0])), axis=1)
    df['point2'] = df.apply(lambda x: Point((x['long'][1], x['lat'][1])), axis=1)
    df['point1_name'] = df.node_id.str[0]
    df['point2_name'] = df.node_id.str[1]

    part1 = df[['point1', 'point1_name']]
    part2 = df[['point2', 'point2_name']]
    part1.columns = ['geometry', 'name']
    part2.columns = ['geometry', 'name']
    points = [part1, part2]
    points = concat_gdf(points)
    points = points.drop_duplicates()
    points.reset_index(drop=True, inplace=True)

    return df, points


def load_bus_regions(onshore_path, offshore_path):
    """Load pypsa-eur on- and offshore regions and concat."""

    bus_regions_offshore = gpd.read_file(offshore_path)
    bus_regions_onshore = gpd.read_file(onshore_path)
    bus_regions = concat_gdf([bus_regions_offshore, bus_regions_onshore])
    bus_regions = bus_regions.dissolve(by='name', aggfunc='sum')
    bus_regions = bus_regions.reset_index()

    return bus_regions


def points2buses(input_points, bus_regions):
    """Map gas network points to network buses depending on bus region."""

    points = input_points.copy()
    points['bus'] = None
    buses_list = set(bus_regions.name)
    for bus in buses_list:
        mask = bus_regions[bus_regions.name == bus]
        index = gpd.clip(points, mask).index
        points.loc[index, 'bus'] = bus

    return points


def build_gas_network_topology(df, points2buses):
    """Create gas network between pypsa buses.

    Parameters
    ----------
    df : pd.DataFrame
        gas network data
    points2buses_map : pd.DataFrame
        mapping of gas network points to pypsa buses

    Returns
    -------
    gas_connections : pd.DataFrame
        gas network connecting pypsa buses
    """

    tmp_df = points2buses[['bus', 'name']]

    tmp_df.columns = ['buses_start', 'name']
    gas_connections = df.merge(tmp_df, left_on='point1_name', right_on='name')
    
    tmp_df.columns = ['buses_destination', 'name']
    gas_connections = gas_connections.merge(tmp_df, left_on='point2_name', right_on='name')
    
    # drop all pipes connecting the same bus
    gas_connections = gas_connections[gas_connections.buses_start != gas_connections.buses_destination]
    gas_connections.reset_index(drop=True, inplace=True)
    gas_connections.drop(['point1', 'point2'], axis=1, inplace=True)

    return gas_connections


def check_missing(nodes, gas_connections):
    """Check which nodes are not connected to the gas network."""

    start_buses = gas_connections.buses_start.dropna().unique()
    end_buses = gas_connections.buses_destination.dropna().unique()

    missing_start = nodes[[bus not in start_buses for bus in nodes]]
    missing_end = nodes[[bus not in end_buses for bus in nodes]]

    logger.info(f"- The following buses are missing in gas network data as a start bus:"
                f"\n {', '.join(map(str, missing_start))} \n"
                f"- The following buses are missing in gas network data as an end bus:"
                f"\n {', '.join(map(str, missing_end))} \n"
                f"- The following buses are missing completely:"
                f"\n {', '.join(map(str, missing_start.intersection(missing_end)))}")


def clean_dataset(nodes, gas_connections):
    """Convert units and save only necessary data."""

    check_missing(nodes, gas_connections)

    determine_pipe_capacity(gas_connections)
    
    cols = [
        'is_bothDirection',
        'capacity_recalculated',
        'buses_start',
        'buses_destination',
        'id',
        'length_km'
    ]
    clean_pipes = gas_connections[cols].dropna()


    # convert GW -> MW
    clean_pipes.loc[:, 'capacity_recalculated'] *= 1e3

    # rename columns
    to_rename = {
        'capacity_recalculated': 'pipe_capacity_MW',
        'buses_start': 'bus0',
        'buses_destination': 'bus1'
    }
    clean_pipes.rename(columns=to_rename, inplace=True)

    return clean_pipes


def diameter2capacity(pipe_diameter_mm):
    """Calculate pipe capacity based on diameter.

    20 inch (500 mm)  50 bar -> 1.5   GW CH4 pipe capacity (LHV)
    24 inch (600 mm)  50 bar -> 5     GW CH4 pipe capacity (LHV)
    36 inch (900 mm)  50 bar -> 11.25 GW CH4 pipe capacity (LHV)
    48 inch (1200 mm) 80 bar -> 21.7  GW CH4 pipe capacity (LHV)

    Based on p.15 of https://gasforclimate2050.eu/wp-content/uploads/2020/07/2020_European-Hydrogen-Backbone_Report.pdf
    """

    # slopes definitions
    m0 = (5 - 1.5) / (600 - 500)
    m1 = (11.25 - 5) / (900 - 600)
    m2 = (21.7 - 11.25) / (1200 - 900)

    # intercept
    a0 = -16
    a1 = -7.5
    a2 = -20.1

    if pipe_diameter_mm < 500:
        return np.nan
    elif pipe_diameter_mm < 600:
        return a0 + m0 * pipe_diameter_mm
    elif pipe_diameter_mm < 900:
        return a1 + m1 * pipe_diameter_mm
    else:
        return a2 + m2 * pipe_diameter_mm


def determine_pipe_capacity(gas_network):
    """Check pipe capacity depending on diameter and pressure."""

    gas_network["capacity_recalculated"] = gas_network.diameter_mm.apply(diameter2capacity)
    
    # if pipe capacity smaller than 1.5 GW take original pipe capacity
    low_cap = gas_network.Capacity_GWh_h < 1.5
    gas_network.loc[low_cap, "capacity_recalculated"] = gas_network.loc[low_cap, "capacity_recalculated"].fillna(gas_network.loc[low_cap, "Capacity_GWh_h"])
    
    # for pipes without diameter assume 500 mm diameter
    gas_network["capacity_recalculated"].fillna(1.5, inplace=True)
    
    # for nord stream take orginal data
    nord_stream = gas_network[gas_network.max_pressure_bar==220].index
    gas_network.loc[nord_stream, "capacity_recalculated"] = gas_network.loc[nord_stream, "Capacity_GWh_h"]


if __name__ == "__main__":

    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_gas_network',
            network='elec', simpl='', clusters='37',
            lv='1.0', opts='', planning_horizons='2020',
            sector_opts='168H-T-H-B-I')

    logging.basicConfig(level=snakemake.config['logging_level'])

    # import gas network data
    gas_network, points = load_gas_network(snakemake.input.gas_network)

    # get clustered bus regions
    bus_regions = load_bus_regions(
        snakemake.input.regions_onshore,
        snakemake.input.regions_offshore
    )
    nodes = pd.Index(bus_regions.name.unique())

    # map gas network points to network buses
    points2buses_map = points2buses(points, bus_regions)

    #  create gas network between pypsa nodes
    gas_connections = build_gas_network_topology(gas_network, points2buses_map)

    gas_connections = clean_dataset(nodes, gas_connections)

    gas_connections.to_csv(snakemake.output.clustered_gas_network)
