#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Builds natural gas network based on data from the SciGRID Gas project
(https://www.gas.scigrid.de/) and ENTSOG capacity map
(https://www.entsog.eu/sites/default/files/2019-10/Capacities%20for%20Transmission%20Capacity%20Map%20RTS008_NS%20-%20DWH_final.xlsx)

Relevant Settings
-----------------
.. code:: yaml
    sector:


Inputs
------
gas network data from SciGRID gas and ENTSOG:
    - IGGINL='data/gas_network/IGGINL_PipeSegments.csv'
    - entsog_2019='data/gas_network/entsog_2019_dataset.csv'
    - EMAP='data/gas_network/EMAP_Raw_PipeSegments.csv'

Outputs
-------
combined gas network data for corresponding PyPSA-Eur-Sec network
    - gas_network='resources/gas_network_{clusters}.csv'

Description
-----------

"""

import logging
logger = logging.getLogger(__name__)
from _helpers import configure_logging

import re
import pandas as pd
import numpy as np
import json

from shapely.geometry import LineString,Point

from sklearn.linear_model import Lasso
from sklearn.metrics import mean_absolute_error

# helper functions ------------------------------------------------------------
def map_entsog_point_to_IGG(entsog_index, IGGINL):
    """
    maps ENTSOG point to closest IGG pipe segment which is connecting the same
    countries
    """
    # missing countries are labelled with "XX" in IGG dataset
    countries = entsog.loc[entsog_index,["From", "To"]].to_list()+["XX"]
    # do not consider direction of pipe
    igg_points = IGGINL[(IGGINL["from"].isin(countries)) & (IGGINL["to"].isin(countries))]

    # get from the IGG pipes connecting the same countries as ENTSOG pipe the closest
    closest = (igg_points.mid_Point.apply(lambda x: x.distance(entsog.Point[entsog_index]))
               .sort_values().iloc[:1].reset_index()
               .rename(columns={"index":"IGG_index", "mid_Point":"distance"}))
    closest["entsog_index"] = entsog_index
    # rename back to original IGG index
    closest["IGG_index"] = closest["IGG_index"].apply(lambda x: x% len(IGGINL))

    return closest


def map_EMAP(series, EMAP_Raw, class_dict={'S': 400, 'M': 700, 'L': 1000}, threshold=0.6,):
    """
    maps EMAP pipe diameter classes to closest gas network pipes with uncertain
    pipe diameter
    if the distance is larger than the threshold distance, original values are
    kept
    """
    if series.loc["uncertain_diameter_mm"]!=0:

        distance = pd.DataFrame()
        distance[0] = EMAP_Raw.Point0.apply(lambda x: x.distance(series.Point0))
        distance[1] = EMAP_Raw.Point1.apply(lambda x: x.distance(series.Point1))
        distance[2] = EMAP_Raw.Point0.apply(lambda x: x.distance(series.mid_Point))
        distance[3] = EMAP_Raw.Point1.apply(lambda x: x.distance(series.mid_Point))

        average_dist = distance.sum(axis=1).sort_values().iloc[:1] / (len(distance.columns))

        if all(average_dist < threshold):
            series['EMAP_class'] = EMAP_Raw.loc[average_dist.index, 'pipe_class_EMap'].values[0]
            series['diameter_mm'] = EMAP_Raw.loc[average_dist.index, 'pipe_class_EMap'].map(class_dict).values[0]
    return series

# main functions --------------------------------------------------------------
def prepare_datasets():
    """
    this function prepares the following 3 dataset to be used in pypsa-eur-sec:

        (1) Scigrid gas data set IGGINL (https://doi.org/10.5281/zenodo.4288440)
        (2) ENTSOG capacity map
        (3) SciGRID gas data set EMAP

    """

    # (1) read and prepocess IGGINL dataset -----------------------------------
    IGGINL = pd.read_csv(snakemake.input.IGGINL, sep=';')
    # convert json format to columns
    IGGINL = pd.concat([IGGINL, IGGINL.param.apply(eval).apply(pd.Series)],
                       axis=1)
    # uncertainty parameters
    uncertainty_parameters = ['max_cap_M_m3_per_d', 'diameter_mm', 'max_pressure_bar']
    # convert json to columns and rename to avoid duplicate column names
    uncertainty_df = (IGGINL.uncertainty.apply(eval).apply(pd.Series)[uncertainty_parameters]
                      .rename(columns=lambda x: "uncertain_" + x))
    IGGINL = pd.concat([IGGINL, uncertainty_df], axis=1)

    # add from to country
    IGGINL['from'] = IGGINL.country_code.apply(lambda x: x.split("'")[1]).str.strip()
    IGGINL['to'] = IGGINL.country_code.apply(lambda x: x.split("'")[3]).str.strip()

    # get the buses
    IGGINL["bus0"] = IGGINL.node_id.apply(lambda x: x.split("'")[1])
    IGGINL["bus1"] = IGGINL.node_id.apply(lambda x: x.split("'")[3])

    # combine long lat to shapely point, take midpoint of pipe segment
    long = IGGINL.long.apply(lambda x: sum(eval(x)) / len(eval(x)))
    lat = IGGINL.lat.apply(lambda x: sum(eval(x)) / len(eval(x)))
    IGGINL['mid_Point'] = (pd.concat([long, lat], axis=1)
                       .apply(lambda x: Point(x['long'], x['lat']), axis=1))
    IGGINL['Point0'] = IGGINL.apply(lambda x: Point(eval(x['long'])[0], eval(x['lat'])[0]), axis=1)
    IGGINL['Point1'] = IGGINL.apply(lambda x: Point(eval(x['long'])[-1], eval(x['lat'])[-1]), axis=1)


    # convert capacity from 1e6*m^3/day -> MWh/h
    # TODO: units NOT really clear in documentation
    # documentation p.18:
    # max_cap_M_m3_per_d: The maximum annual gas volume that the pipe can transmit in units of [Mm 3 d −1 ].
    # documentation p.50
    # daily gas flow capacity

    # varibales:
    # v: velocity, Q: volumetric_flow, d: pipe diameter, A: cross-sectional area
    # v = Q /A =  Q / (pi * (d/2)^2)
    # to get sensefull gas flow velocities (v=5-10 m/s, sometimes 20m/s) volumetric flow should be annual
    velocity = IGGINL.max_cap_M_m3_per_d * 1e6 / 8760 / 60 / 60 / (np.pi * (IGGINL.diameter_mm * 1e-3 * 0.5)**2)

    # specific gas constant methane R_s=518.4 J/(kgK)
    R_s = 518.4
    # temperature [Kelvin] (assuming 10°Celsius)
    T = 10 + 273.15
    # density [kg/m^3]= pressure [kg/ms^2] / (T * R_s), 1 bar = 1e5 kg/(ms^2)
    density = IGGINL.max_pressure_bar * 1e5 / (T * R_s)
    # mass flow [kg/ h], Mega = 1e6,
    mass_flow = IGGINL.max_cap_M_m3_per_d * 1e6 / 8760  * density
    # gross calorific value (GCV in ENTSOT table) [kWh/kg]
    gcv_lgas = 38.3 / 3.6
    gcv_hgas = 47.3 / 3.6
    # energy cap [MW] = mass_flow [kg/h] * gcv [kWh/kg] * 1e-3
    energy_cap = mass_flow * 1e-3
    energy_cap.loc[IGGINL.is_H_gas==1] *= gcv_hgas
    energy_cap.loc[IGGINL.is_H_gas!=1] *= gcv_lgas
    IGGINL['max_capacity'] = energy_cap


    # (2) read and preprocess ENTSOG data -------------------------------------
    entsog = pd.read_csv(snakemake.input.entsog_2019)
    entsog.drop(['From_ID', 'To_ID'], axis=1, inplace=True)
    # group parallel pipes and take maximum pipe capacity
    entsog_wrapping = entsog.groupby(['long', 'lat', 'From', 'To']).max()['Capacity'].reset_index()
    # add shapely object
    entsog_wrapping['Point'] = entsog_wrapping.apply(lambda x: Point(x['long'], x['lat']), axis=1)
    # convert GWh/day to MW
    entsog_wrapping["Capacity"] *= 1e3

    # (3) read and preprocess EMAP data ---------------------------------------
    EMAP_Raw = pd.read_csv(snakemake.input.EMAP, sep=';')
    # convert json format to columns
    EMAP_Raw = pd.concat([EMAP_Raw, EMAP_Raw.param.apply(eval).apply(pd.Series)],
                       axis=1)
    # fill missing pipe size (["S", "M", "L"]) values with "A"
    EMAP_Raw.pipe_class_EMap = EMAP_Raw.pipe_class_EMap.fillna('A')
    # add shapely object
    EMAP_Raw['Point0'] = EMAP_Raw.apply(lambda x: Point(eval(x['long'])[0], eval(x['lat'])[0]), axis=1)
    EMAP_Raw['Point1'] = EMAP_Raw.apply(lambda x: Point(eval(x['long'])[-1], eval(x['lat'])[-1]), axis=1)

    return IGGINL, entsog_wrapping, EMAP_Raw


def train_lasso(IGG, alpha=0.001):
    """
    trains lasso regression with unapproximated data of IGG with known
    pipe diameter, pressure and capacity

    normal lasso method is choosen
    """
    # ------------preprocessing----------------
    # find all pipe that have not approximated diameter, pressure and capacity data
    all_data_i = (IGG.loc[:,IGG.columns.str.contains("uncertain_")]==0).all(axis=1)
    train_data = IGG[all_data_i].reset_index(drop=True)

    # capacity depends on squared diameter -> add diameter^2 to training data
    train_data['diameter_squared'] = train_data.diameter_mm ** 2

    # -------------start training--------------
    logger.info('training lasso')
    rg_model_normal = Lasso(alpha=alpha)
    rg_model_normal.fit(train_data.diameter_mm.values.reshape(-1, 1),
                        train_data.max_cap_M_m3_per_d)
    train_data['predict_normal'] = rg_model_normal.predict(
        train_data.diameter_mm.values.reshape(-1, 1))

    # calculate mean absolute error (MAE)
    MAE = str(round(mean_absolute_error(train_data.max_cap_M_m3_per_d,
                                               train_data.predict_normal), 3))

    logger.info('sucessful training lasso regression, mean absoulte error (MAE) = {}' \
                .format(MAE))


    return rg_model_normal


def add_entsog_capacity(IGGINL, entsog, distance_threshold=0.5):
    '''
    merges IGGINL and entsog crossborder pipe capacities, currently pipe
    directions are not considered

    '''
    gas_network = IGGINL.copy()

    # find for every entsog point closest IGG point
    distance = (pd.concat([map_entsog_point_to_IGG(i, IGGINL) for i in entsog.index])
               .groupby("IGG_index").min())

    # get all points within threshold
    distance = distance[distance.distance<distance_threshold].reset_index()
    # set capacitiy
    IGGINL.loc[distance.IGG_index, 'max_capacity'] = entsog.loc[distance.entsog_index, "Capacity"]
    IGGINL['distance_to_capacity_point'] = np.nan
    IGGINL.loc[distance.IGG_index, 'distance_to_capacity_point'] = distance.set_index("IGG_index")["distance"]

    logger.info('adding  {} pipe capacities from ENTSOG  '
                    .format(len(distance)))

    return gas_network


def add_EMAP_diameter(gas_network, EMAP, threshold=0.6):
    """
    add EMAP diameter to the combined data set gas_network for diameters with
    uncertainty
    """
    # calculate mean value of each class with original diameter data
    gas_network_diameter = gas_network.loc[gas_network.uncertain_diameter_mm==0]["diameter_mm"]
    # get mean IGG diameter for pipe classes S,M,L
    gas_network_mean_diameter_s = gas_network_diameter[gas_network_diameter < 600].mean()
    gas_network_mean_diameter_m = gas_network_diameter[(gas_network_diameter >= 600) & (gas_network_diameter < 900)].mean()
    gas_network_mean_diameter_l = gas_network_diameter[gas_network_diameter >= 900].mean()
    pipe_diameter_dict = {'S': gas_network_mean_diameter_s,
                          'M': gas_network_mean_diameter_m,
                          'L': gas_network_mean_diameter_l}

    # filter on EMAP, length>50, only keep S M L
    EMAP = EMAP[EMAP.length_km > 50]
    EMAP = EMAP[EMAP.pipe_class_EMap.isin(['S', 'M', 'L'])]
    EMAP = EMAP.reset_index(drop=True)

    # start matching
    gas_network['EMAP_class'] = np.nan
    gas_network = gas_network.apply(lambda x: map_EMAP(x, EMAP,
                                                       pipe_diameter_dict,
                                                       threshold=threshold),
                                    axis=1)

    logger.info('adding  {} pipe diameters from EMAP  '
                .format(gas_network["EMAP_class"].notna().sum()))

    return gas_network


def filling_with_lasso(gas_network, regression_model):
    """
    fills uncertain values with own lasso regression model

    if diameter of a pipe is still missing, use diameter data from
    diameter_mm of the pipe
    """

    uncertain = gas_network['uncertain_max_cap_M_m3_per_d']!=0
    minimum_value = gas_network[~uncertain]['max_cap_M_m3_per_d'].min()

    gas_network.capacity_nan = gas_network.apply(
        lambda x: regression_model.predict(np.array([x['diameter_nan']]).reshape(-1, 1))[0]
        if (np.isnan(x['capacity_nan'])) & (not np.isnan(x['diameter_nan'])) else x['capacity_nan'], axis=1)

    #remove extremely small value
    gas_network.capacity_nan = gas_network.capacity_nan.apply(lambda x: np.nan if x < minimum_value else x)
    logger.info('finish filling with lasso')
    return gas_network
#%%
if __name__ == "__main__":

    # for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake
        snakemake = MockSnakemake(
            wildcards=dict(network='elec', simpl='', clusters='37', lv='1.0',
                           opts='', planning_horizons='2020',
                           sector_opts='168H-T-H-B-I'),

            input=dict(IGGINL='data/gas_network/IGGINL_PipeSegments.csv',
                       entsog_2019='data/gas_network/entsog_2019_dataset.csv',
                       EMAP='data/gas_network/EMAP_Raw_PipeSegments.csv'
                       ),
            output=dict(gas_network='resources/gas_network_{clusters}.csv'),
        )
        import yaml
        with open('config.yaml', encoding='utf8') as f:
            snakemake.config = yaml.safe_load(f)


    # prepare the data sets
    IGGINL, entsog, EMAP = prepare_datasets()

    # train lasso regression
    regression_model = train_lasso(IGGINL)

    # add crossborder capacities from ENTSOG
    gas_network = add_entsog_capacity(IGGINL, entsog)

    # TODO ------------------------------------------------------
    # add pipe diameters from EMAP
    gas_network = add_EMAP_diameter(gas_network, EMAP)

    # fill other missing values with lasso regression model
    IGGINL = filling_with_lasso(IGGINL, regression_model)

    # IGGINL = node_capacity_spread(IGGINL)
    # clean_save(IGGINL, output_path)
