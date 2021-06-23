#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Builds clustered natural gas network based on data from:
    [1] the SciGRID Gas project
        (https://www.gas.scigrid.de/)
    [2] ENTSOG capacity map
        (https://www.entsog.eu/sites/default/files/2019-10/Capacities%20for%20Transmission%20Capacity%20Map%20RTS008_NS%20-%20DWH_final.xlsx)

Relevant Settings
-----------------
.. code:: yaml
    sector:


Inputs
------
gas network data from SciGRID gas and ENTSOG:
    - gas_network="data/gas_network/gas_network_dataset.csv",
      combined gas network data set from [1] and [2]
    - country_shapes=pypsaeur("resources/country_shapes.geojson"),
    - regions_onshore=pypsaeur("resources/regions_onshore_elec_s{simpl}_{clusters}.geojson"),
    - regions_offshore=pypsaeur("resources/regions_offshore_elec_s{simpl}_{clusters}.geojson")

Outputs
-------
clustered gas network data for corresponding PyPSA-Eur-Sec network
    - gas_network='resources/gas_network_{clusters}.csv'

Description
-----------

"""
import geoplot
import geoplot.crs as gcrs
import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__name__)

import re
import pandas as pd
import json

from shapely.geometry import LineString,Point

import geopandas as gpd
import numpy as np

#-----------------##########################################################
# helper functions #
#-----------------#
def concat_gdf(gdf_list, crs = 'EPSG:4326'):
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


#-----------------############################################################
#  main functions #
#-----------------#
def preprocessing(df_path):
    """Load  and format gas network data."""
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


def load_region(onshore_path, offshore_path):
    """Load pypsa-eur on- and offshore regions and concat."""
    buses_region_offshore = gpd.read_file(offshore_path)
    buses_region_onshore = gpd.read_file(onshore_path)
    buses_region = concat_gdf([buses_region_offshore, buses_region_onshore])
    buses_region = buses_region.dissolve(by='name', aggfunc='sum')
    buses_region = buses_region.reset_index()

    return buses_region


def create_points2buses_map(input_points, buses_region):
    """Map gas network points to network buses depending on bus region."""
    points = input_points.copy()
    points['bus'] = None
    buses_list = set(buses_region.name)
    for bus in buses_list:
        mask = buses_region[buses_region.name == bus]
        index = gpd.clip(points, mask).index
        points.loc[index, 'bus'] = bus

    return points


def create_cross_regions_network(df, points2buses_map):
    """Create gas network between pypsa buses.

    Input:
        df               :  gas network data (pd.DataFrame)
        points2buses_map : map gas network points to pypsa buses (pd.DataFrame)

    Return:
        cross_buses_gas_network : gas network connecting pypsa buses
                                  (pd.DataFrame)
    """
    tmp_df = points2buses_map[['bus', 'name']]
    tmp_df.columns = ['buses_start','name']
    cross_buses_gas_network = df.merge(tmp_df, left_on='point1_name',
                                       right_on='name')
    tmp_df.columns = ['buses_destination', 'name']
    cross_buses_gas_network = cross_buses_gas_network.merge(tmp_df,
                                                            left_on='point2_name',
                                                            right_on='name')
    # drop all pipes connecting the same bus
    cross_buses_gas_network = cross_buses_gas_network[cross_buses_gas_network.buses_start \
                                                      != cross_buses_gas_network.buses_destination]
    cross_buses_gas_network.reset_index(drop=True, inplace=True)
    cross_buses_gas_network.drop(['point1','point2'], axis=1, inplace=True)

    return cross_buses_gas_network


def check_missing(nodes, cross_buses_gas_network):
    """Check which nodes are not connected to the gas network."""
    missing0 = nodes[[bus not in cross_buses_gas_network.buses_start.dropna().unique()
                     for bus in nodes]]
    missing1 = nodes[[bus not in cross_buses_gas_network.buses_destination.dropna().unique()
                 for bus in nodes]]
    logger.info("\n - The following buses are missing in gas network data as a start bus: \n {} \n"
                "- The following buses are missing in gas network data as an end bus: \n {} \n "
                "- The following buses are missing completely: \n {}"
                .format(', '.join(map(str, missing0)),
                        ', '.join(map(str, missing1)),
                        ', '.join(map(str, missing0.intersection(missing1)))))


def clean_dataset(cross_buses_gas_network):
    """Convert units and save only necessary data."""
    inspect_pipe_capacity(cross_buses_gas_network)
    cols = ['is_bothDirection', 'capacity_recalculated','buses_start',
            'buses_destination', 'id', 'length_km']
    clean_pipes = cross_buses_gas_network[cols].dropna()


    # convert GW -> MW
    clean_pipes.loc[:, 'capacity_recalculated'] *= 1e3
    # rename columns
    clean_pipes.rename(columns={'capacity_recalculated': 'pipe_capacity_MW',
                        'buses_start': 'bus0',
                        'buses_destination': 'bus1'}, inplace=True)
    return clean_pipes


def recalculate_pipe_capacity(pipe_diameter_mm):
    """Calculate pipe capacity based on diameter.

    20 inch (500 mm)  50 bar -> 1.5   GW CH4 pipe capacity (LHV)
    24 inch (600 mm)  50 bar -> 5     GW CH4 pipe capacity (LHV)
    36 inch (900 mm)  50 bar -> 11.25 GW CH4 pipe capacity (LHV)
    48 inch (1200 mm) 80 bar -> 21.7  GW CH4 pipe capacity (LHV)

    Based on p.15 of (https://gasforclimate2050.eu/wp-content/uploads/2020/07/2020_European-Hydrogen-Backbone_Report.pdf"""
    # slope
    m0 = (5-1.5) / (600-500)
    m1 = (11.25-5)/(900-600)
    m2 = (21.7-11.25)/(1200-900)

    if pipe_diameter_mm<500:
        return np.nan
    if pipe_diameter_mm<600 and pipe_diameter_mm>=500:
        return -16 + m0 * pipe_diameter_mm
    if pipe_diameter_mm<900 and pipe_diameter_mm>=600:
        return -7.5 + m1 * pipe_diameter_mm
    else:
        return -20.1 + m2 * pipe_diameter_mm

def inspect_pipe_capacity(gas_network):
    """Check pipe capacity depending on diameter and pressure."""
    gas_network["capacity_recalculated"] = gas_network.diameter_mm.apply(recalculate_pipe_capacity)
    low_cap = gas_network.Capacity_GWh_h < 1.5
    # if pipe capacity smaller than 1.5 GW take original pipe capacity
    gas_network.loc[low_cap, "capacity_recalculated"] = gas_network.loc[low_cap, "capacity_recalculated"].fillna(gas_network.loc[low_cap,"Capacity_GWh_h"])
    gas_network["capacity_recalculated"].fillna(1.5, inplace=True)
    # nord stream take orginal data
    nord_stream = gas_network[gas_network.max_pressure_bar==220].index
    gas_network.loc[nord_stream, "capacity_recalculated"] = gas_network.loc[nord_stream, "Capacity_GWh_h"]

# ----------- VISULAISATION --------------------------------------------------
def create_view_object(cbgn_no_duplicate,buses_region):
    """Create object to view gas network data."""
    cbgn_no_duplicate=cbgn_no_duplicate.merge(buses_region,left_on='buses_start',right_on='name')
    cbgn_no_duplicate=cbgn_no_duplicate.merge(buses_region,left_on='buses_destination',right_on='name')

    cbgn_no_duplicate.geometry_x=cbgn_no_duplicate.geometry_x.apply(lambda x: x.centroid)
    cbgn_no_duplicate.geometry_y=cbgn_no_duplicate.geometry_y.apply(lambda x: x.centroid)
    cbgn_no_duplicate['geometry']=list(zip(cbgn_no_duplicate['geometry_x'],
                                           cbgn_no_duplicate['geometry_y']))

    final = cbgn_no_duplicate[['buses_start', 'buses_destination',
                               'Capacity_GWh_h', 'geometry']]
    final['geometry'] = final['geometry'].apply(LineString)
    final=gpd.GeoDataFrame(final,crs='EPSG:4326')

    return final


def view(cbgn_no_duplicate, buses_region, shapes_path):
    """Plot gas network."""
    final = create_view_object(cbgn_no_duplicate,buses_region)

    eu=gpd.read_file(shapes_path)
    ax = geoplot.webmap(eu, projection=gcrs.WebMercator(), figsize=(20,20),
                        alpha=0.5)
    geoplot.choropleth(buses_region, hue='name',ax=ax, alpha=0.2,
                       edgecolor='red', linewidth=2)
    geoplot.sankey( final, scale='Capacity_GWh_h', hue='Capacity_GWh_h',
                   cmap='viridis', ax=ax, legend=True, legend_var='hue')
    plt.savefig("../graphics/clustered-gas-network_{}.pdf".format(snakemake.wildcards.clusters),
                bbox_inches='tight', pad_inches=0.1)

#%%
if __name__ == "__main__":

    # for testing
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('build_gas_network',
                                   network='elec', simpl='', clusters='37',
                                   lv='1.0', opts='', planning_horizons='2020',
                                   sector_opts='168H-T-H-B-I')

    logging.basicConfig(level=snakemake.config['logging_level'])

    # import gas network data
    gas_network, points = preprocessing(snakemake.input.gas_network)

    # get clustered bus regions
    buses_region = load_region(snakemake.input.regions_onshore,
                               snakemake.input.regions_offshore)
    nodes = pd.Index(buses_region.name.unique())

    # map gas network points to network buses
    points2buses_map = create_points2buses_map(points, buses_region)
    # create gas network between pypsa nodes
    cross_buses_gas_network = create_cross_regions_network(gas_network,
                                                           points2buses_map)

    # view(cross_buses_gas_network, buses_region, snakemake.input.country_shapes)

    # check which buses are not connected in gas network
    check_missing(nodes, cross_buses_gas_network)

    # convert units and save only needed data
    gas_pipes = clean_dataset(cross_buses_gas_network)
    gas_pipes.to_csv(snakemake.output.clustered_gas_network)
