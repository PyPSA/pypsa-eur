#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 14:57:21 2020

@author: bw0928

*****************************************************************************
This script calculates cost-energy_saving-curves for retrofitting
for the EU-37 countries, based on the building stock data from hotmaps and
the EU building stock database
*****************************************************************************

Structure:

    (1) set assumptions and parameters
    (2) read and prepare data
    (3) calculate (€-dE-curves)
    (4) save in csv

*****************************************************************************
"""

import pandas as pd
import matplotlib.pyplot as plt

#%% ************ FUCNTIONS ***************************************************

# windows ---------------------------------------------------------------
def window_limit(l, window_assumptions):
    """
    define limit u value from which on window is retrofitted
    """
    m = (window_assumptions.diff()["u_limit"] /
         window_assumptions.diff()["strength"]).dropna().iloc[0]
    a = window_assumptions["u_limit"][0] - m * window_assumptions["strength"][0]
    return m*l + a

def u_retro_window(l, window_assumptions):
    """
    define retrofitting value depending on renovation strength
    """
    m = (window_assumptions.diff()["u_value"] /
         window_assumptions.diff()["strength"]).dropna().iloc[0]
    a = window_assumptions["u_value"][0] - m * window_assumptions["strength"][0]
    return max(m*l + a, 0.8)

def window_cost(u, cost_retro, window_assumptions):
    """
    get costs for new windows depending on u value

    """
    m = (window_assumptions.diff()["cost"] /
         window_assumptions.diff()["u_value"]).dropna().iloc[0]
    a = window_assumptions["cost"][0] - m * window_assumptions["u_value"][0]
    window_cost = m*u + a
    if annualise_cost:
        window_cost = window_cost * interest_rate / (1 - (1 + interest_rate)
                      ** -cost_retro.loc["Windows", "life_time"])
    return window_cost

#  functions for intermediate steps (~l, ~area)  -----------------------------
def calculate_new_u(u_values, l, l_weight, k=0.035):
    """
    calculate U-values after building retrofitting, depending on the old
    U-values (u_values).
    They depend for the components Roof, Wall, Floor on the additional
    insulation thickness (l), and the weighting for the corresponding
    component (l_weight).
    Windows are renovated to new ones with U-value (function: u_retro_window(l))
    only if the are worse insulated than a certain limit value
    (function: window_limit).

    Parameters
    ----------
    u_values: pd.DataFrame
    l: string
    l_weight: pd.DataFrame (component, weight)
    k: thermal conductivity

    """
    return u_values.apply(lambda x:
                                 k / ((k / x.value) +
                                      (float(l) * l_weight.loc[x.type][0]))
                                 if x.type!="Windows"
                             else (min(x.value, u_retro_window(float(l), window_assumptions))
                                   if x.value>window_limit(float(l), window_assumptions) else x.value),
                             axis=1)

def calculate_dE(u_values, l, average_surface_w):
    """
    returns energy demand after retrofit (per unit of unrefurbished energy
    demand) depending on current and retrofitted U-values, this energy demand
    is weighted depending on the average surface of each component for the
    building type of the assumend subsector
    """
    return u_values.apply(lambda x: x[l] / x.value *
                                     average_surface_w.loc[x.assumed_subsector,
                                                           x.type],
                                     axis=1)


def calculate_costs(u_values, l, cost_retro, average_surface):
    """
    returns costs for a given retrofitting strength weighted by the average
    surface/volume ratio of the component for each building type
    """
    return u_values.apply(lambda x: (cost_retro.loc[x.type, "cost_var"] *
                                     100 * float(l) * l_weight.loc[x.type][0]
                                     + cost_retro.loc[x.type, "cost_fix"]) *
                              average_surface.loc[x.assumed_subsector, x.type] /
                              average_surface.loc[x.assumed_subsector, "surface"]
                              if x.type!="Windows"
                              else (window_cost(x[l], cost_retro, window_assumptions) *
                                               average_surface.loc[x.assumed_subsector, x.type] /
                                               average_surface.loc[x.assumed_subsector, "surface"]
                                    if x.value>window_limit(float(l), window_assumptions) else 0),
                             axis=1)


# ---------------------------------------------------------------------------
def prepare_building_stock_data():
    """
    reads building stock data and cleans up the format, returns
    --------
    u_values:          pd.DataFrame current U-values
    average_surface:   pd.DataFrame (index= building type,
                                     columns = [surface [m],height [m],
                                                components area [m^2]])
    average_surface_w: pd.DataFrame weighted share of the components per
                       building type
    area_tot:          heated floor area per country and sector [Mm²]
    area:              heated floor area [Mm²] for country, sector, building
                       type and period

    """

    building_data = pd.read_csv(snakemake.input.building_stock,
                                usecols=list(range(13)))

    # standardize data
    building_data["type"].replace(
                {'Covered area: heated  [Mm²]': 'Heated area [Mm²]',
                 'Windows ': 'Windows',
                 'Walls ': 'Walls',
                 'Roof ': 'Roof',
                 'Floor ': 'Floor'}, inplace=True)

    building_data.country_code = building_data.country_code.str.upper()
    building_data["subsector"].replace({'Hotels and Restaurants':
                                        'Hotels and restaurants'}, inplace=True)
    building_data["sector"].replace({'Residential sector': 'residential',
                                     'Service sector': 'services'},
                                    inplace=True)
    # extract u-values
    u_values = building_data[(building_data.feature.str.contains("U-values"))
                             & (building_data.subsector != "Total")]

    components = list(u_values.type.unique())

    country_iso_dic = building_data.set_index("country")["country_code"].to_dict()

    # add missing /rename countries
    country_iso_dic.update({'Norway': 'NO',
                            'Iceland': 'IS',
                            'Montenegro': 'ME',
                            'Serbia': 'RS',
                            'Albania': 'AL',
                            'United Kingdom': 'GB',
                            'Bosnia and Herzegovina': 'BA',
                            'Switzerland': 'CH'})

    # heated floor area ----------------------------------------------------------
    area = building_data[(building_data.type == 'Heated area [Mm²]') &
                         (building_data.subsector != "Total")]
    area_tot = area.groupby(["country", "sector"]).sum()
    area = pd.concat([area, area.apply(lambda x: x.value /
                                          area_tot.value.loc[(x.country, x.sector)],
                                          axis=1).rename("weight")],axis=1)
    area = area.groupby(['country', 'sector', 'subsector', 'bage']).sum()
    area_tot.rename(index=country_iso_dic, inplace=True)

    # add for some missing countries floor area from other data sources
    area_missing = pd.read_csv(snakemake.input.floor_area_missing,
                               index_col=[0, 1], usecols=[0, 1, 2, 3],
                               encoding='ISO-8859-1')
    area_tot = area_tot.append(area_missing.unstack(level=-1).dropna().stack())
    area_tot = area_tot.loc[~area_tot.index.duplicated(keep='last')]

    # for still missing countries calculate floor area by population size
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    pop_layout["ct"] = pop_layout.index.str[:2]
    ct_total = pop_layout.total.groupby(pop_layout["ct"]).sum()

    area_per_pop = area_tot.unstack().reindex(index=ct_total.index).apply(lambda x: x / ct_total[x.index])
    missing_area_ct = ct_total.index.difference(area_tot.index.levels[0])
    for ct in missing_area_ct.intersection(ct_total.index):
        averaged_data = pd.DataFrame(
            area_per_pop.value.reindex(map_for_missings[ct]).mean()
            * ct_total[ct],
            columns=["value"])
        index = pd.MultiIndex.from_product([[ct], averaged_data.index.to_list()])
        averaged_data.index = index
        averaged_data["estimated"] = 1
        if ct not in area_tot.index.levels[0]:
            area_tot = area_tot.append(averaged_data, sort=True)
        else:
            area_tot.loc[averaged_data.index] = averaged_data

    # u_values for Poland are missing -> take them from eurostat -----------
    u_values_PL = pd.read_csv(snakemake.input.u_values_PL)
    area_PL = area.loc["Poland"].reset_index()
    data_PL = pd.DataFrame(columns=u_values.columns, index=area_PL.index)
    data_PL["country"] = "Poland"
    data_PL["country_code"] = "PL"
    # data from area
    for col in ["sector", "subsector", "bage"]:
        data_PL[col] = area_PL[col]
    data_PL["btype"] = area_PL["subsector"]

    data_PL_final = pd.DataFrame()
    for component in components:
        data_PL["type"] = component
        data_PL["value"] = data_PL.apply(lambda x: u_values_PL[(u_values_PL.component==component)
                                                               & (u_values_PL.sector==x["sector"])]
                                         [x["bage"]].iloc[0], axis=1)
        data_PL_final = data_PL_final.append(data_PL)

    u_values = pd.concat([u_values,
                          data_PL_final]).reset_index(drop=True)

    # clean data ---------------------------------------------------------------
    # smallest possible today u values for windows 0.8 (passive house standard)
    # maybe the u values for the glass and not the whole window including frame
    # for those types assumed in the dataset
    u_values.loc[(u_values.type=="Windows") & (u_values.value<0.8), "value"] = 0.8
    # drop unnecessary columns
    u_values.drop(['topic', 'feature','detail', 'estimated','unit'],
                  axis=1, inplace=True, errors="ignore")
    #  only take in config.yaml specified countries into account
    countries = ct_total.index
    area_tot = area_tot.loc[countries]

    # average component surface --------------------------------------------------
    average_surface = (pd.read_csv(snakemake.input.average_surface,
                                   nrows=3,
                                   header=1,
                                   index_col=0).rename(
                                       {'Single/two family house': 'Single family- Terraced houses',
                                        'Large apartment house': 'Multifamily houses',
                                        'Apartment house': 'Appartment blocks'},
                                    axis="index")).iloc[:, :6]
    average_surface.columns = ["surface", "height", "Roof",
                               "Walls", "Floor", "Windows"]
    # get area share of component
    average_surface_w = average_surface[components].apply(lambda x: x / x.sum(),
                                                          axis=1)

    return (u_values, average_surface,
            average_surface_w, area_tot, area, country_iso_dic, countries)


def prepare_cost_retro():
    """
    read and prepare retro costs, annualises them if annualise_cost=True
    """
    cost_retro = pd.read_csv(snakemake.input.cost_germany,
                         nrows=4, index_col=0, usecols=[0, 1, 2, 3])
    cost_retro.index = cost_retro.index.str.capitalize()
    cost_retro.rename(index={"Window": "Windows", "Wall": "Walls"}, inplace=True)

    window_assumptions = pd.read_csv(snakemake.input.window_assumptions,
                                     skiprows=[1], usecols=[0,1,2,3], nrows=2)

    if annualise_cost:
        cost_retro[["cost_fix", "cost_var"]] = (cost_retro[["cost_fix", "cost_var"]]
                                                .apply(lambda x: x * interest_rate /
                                                       (1 - (1 + interest_rate)
                                                        ** -cost_retro.loc[x.index,
                                                                           "life_time"])))

    return cost_retro, window_assumptions


def calculate_cost_energy_curve(u_values, l_strength, l_weight, average_surface_w,
                                average_surface, area, country_iso_dic,
                                countries):
    """
    returns energy demand per unit of unrefurbished (dE) and cost for given
    renovation strength (l_strength), data for missing countries is
    approximated by countries with similar building stock (dict:map_for_missings)

    parameter
    -------- input -----------
    u_values:          pd.DataFrame current U-values
    l_strength:        list of strings (strength of retrofitting)
    l_weight:          pd.DataFrame (component, weight)
    average_surface:   pd.DataFrame (index= building type,
                                     columns = [surface [m],height [m],
                                                components area [m^2]])
    average_surface_w: pd.DataFrame weighted share of the components per
                       building type
    area:              heated floor area [Mm²] for country, sector, building
                       type and period
    country_iso_dic:   dict (maps country name to 2-letter-iso-code)
    countries:         pd.Index (specified countries in config.yaml)
    -------- output ----------
    res:               pd.DataFrame(index=pd.MultiIndex([country, sector]),
                                    columns=pd.MultiIndex([(dE/cost), l_strength]))
    """

    energy_saved = u_values[['country', 'sector', 'subsector', 'bage', 'type']]
    costs = u_values[['country', 'sector', 'subsector', 'bage', 'type']]

    for l in l_strength:
        u_values[l] = calculate_new_u(u_values, l, l_weight)
        energy_saved = pd.concat([energy_saved,
                                  calculate_dE(u_values, l, average_surface_w).rename(l)],
                                 axis=1)
        costs = pd.concat([costs,
                           calculate_costs(u_values, l, cost_retro, average_surface).rename(l)],
                          axis=1)

    # energy and costs per country, sector, subsector and year
    e_tot = energy_saved.groupby(['country', 'sector', 'subsector', 'bage']).sum()
    cost_tot = costs.groupby(['country', 'sector', 'subsector', 'bage']).sum()

    # weighting by area -> energy and costs per country and sector
    # in case of missing data first concat
    energy_saved = pd.concat([e_tot, area.weight], axis=1)
    cost_res = pd.concat([cost_tot, area.weight], axis=1)
    energy_saved = (energy_saved.apply(lambda x: x * x.weight, axis=1)
                    .groupby(level=[0, 1]).sum())
    cost_res = (cost_res.apply(lambda x: x * x.weight, axis=1)
                .groupby(level=[0, 1]).sum())

    res = pd.concat([energy_saved[l_strength], cost_res[l_strength]],
                    axis=1, keys=["dE", "cost"])
    res.rename(index=country_iso_dic, inplace=True)

    res = res.reindex(index=countries, level=0)
    # reset index because otherwise not considered countries still in index.levels[0]
    res = res.reset_index().set_index(["country", "sector"])

    # map missing countries
    for ct in pd.Index(map_for_missings.keys()).intersection(countries):
        averaged_data = res.reindex(index=map_for_missings[ct], level=0).mean(level=1)
        index = pd.MultiIndex.from_product([[ct], averaged_data.index.to_list()])
        averaged_data.index = index
        if ct not in res.index.levels[0]:
            res = res.append(averaged_data)
        else:
            res.loc[averaged_data.index] = averaged_data

    return res


# %% **************** MAIN ************************************************
if __name__ == "__main__":
    #  for testing
    if 'snakemake' not in globals():
        import yaml
        import os
        from vresutils.snakemake import MockSnakemake
        snakemake = MockSnakemake(
            wildcards=dict(
                network='elec',
                simpl='',
                clusters='37',
                lv='1',
                opts='Co2L-3H',
                sector_opts="[Co2L0p0-168H-T-H-B-I]"),
            input=dict(
                building_stock="data/retro/data_building_stock.csv",
                u_values_PL="data/retro/u_values_poland.csv",
                tax_w="data/retro/electricity_taxes_eu.csv",
                construction_index="data/retro/comparative_level_investment.csv",
                average_surface="data/retro/average_surface_components.csv",
                floor_area_missing="data/retro/floor_area_missing.csv",
                clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
                cost_germany="data/retro/retro_cost_germany.csv",
                window_assumptions="data/retro/window_assumptions.csv"),
            output=dict(
                retro_cost="resources/retro_cost_elec_s{simpl}_{clusters}.csv",
                floor_area="resources/floor_area_elec_s{simpl}_{clusters}.csv")
        )
        with open('config.yaml', encoding='utf8') as f:
            snakemake.config = yaml.safe_load(f)

#  ******** (1) ASSUMPTIONS - PARAMETERS **********************************
    retro_opts =  snakemake.config["sector"]["retrofitting"]
    interest_rate = retro_opts["interest_rate"]
    annualise_cost = retro_opts["annualise_cost"]  # annualise the investment costs
    tax_weighting = retro_opts["tax_weighting"]   # weight costs depending on taxes in countries
    construction_index = retro_opts["construction_index"]   # weight costs depending on labour/material costs per ct
    # additional insulation thickness, determines maximum possible savings
    l_strength = retro_opts["l_strength"]

    k = 0.035   # thermal conductivity standard value
    # strenght of relative retrofitting depending on the component
    # determined by historical data of insulation thickness for retrofitting
    l_weight = pd.DataFrame({"weight": [1.95, 1.48, 1.]},
                            index=["Roof", "Walls", "Floor"])


    # mapping missing countries by neighbours
    map_for_missings = {
        "AL": ["BG", "RO", "GR"],
        "BA": ["HR"],
        "RS": ["BG", "RO", "HR", "HU"],
        "MK": ["BG", "GR"],
        "ME": ["BA", "AL", "RS", "HR"],
        "CH": ["SE", "DE"],
        "NO": ["SE"],
        }

# %% ************ (2) DATA ***************************************************

    # building stock data -----------------------------------------------------
    (u_values, average_surface, average_surface_w,
     area_tot, area, country_iso_dic, countries) = prepare_building_stock_data()

    # costs for retrofitting -------------------------------------------------
    cost_retro, window_assumptions = prepare_cost_retro()

    # weightings of costs
    if construction_index:
        cost_w = pd.read_csv(snakemake.input.construction_index,
                             skiprows=3, nrows=32, index_col=0)
        # since German retrofitting costs are assumed
        cost_w = ((cost_w["2018"] / cost_w.loc["Germany", "2018"])
                  .rename(index=country_iso_dic))

    if tax_weighting:
        tax_w = pd.read_csv(snakemake.input.tax_w,
                            header=12, nrows=39, index_col=0, usecols=[0, 4])
        tax_w.rename(index=country_iso_dic, inplace=True)
        tax_w = tax_w.apply(pd.to_numeric, errors='coerce').iloc[:, 0]
        tax_w.dropna(inplace=True)

# %% ********** (3) CALCULATE COST-ENERGY-CURVES ****************************

    # for missing weighting of surfaces of building types assume MultiFamily houses
    u_values["assumed_subsector"] = u_values.subsector
    u_values.loc[~u_values.subsector.isin(average_surface.index),
                 "assumed_subsector"] = 'Multifamily houses'

    dE_and_cost =  calculate_cost_energy_curve(u_values, l_strength, l_weight,
                                          average_surface_w, average_surface, area,
                                          country_iso_dic, countries)
    # reset index because otherwise not considered countries still in index.levels[0]
    dE_and_cost = dE_and_cost.reset_index().set_index(["country", "sector"])

    # weights costs after construction index
    if construction_index:
        for ct in list(map_for_missings.keys() - cost_w.index):
            cost_w.loc[ct] = cost_w.reindex(index=map_for_missings[ct]).mean()
        dE_and_cost.cost = dE_and_cost.cost.apply(lambda x: x * cost_w[x.index.levels[0]])

    # weights cost depending on country taxes
    if tax_weighting:
        for ct in list(map_for_missings.keys() - tax_w.index):
            tax_w[ct] = tax_w.reindex(index=map_for_missings[ct]).mean()
        dE_and_cost.cost = dE_and_cost.cost.apply(lambda x: x * tax_w[x.index.levels[0]])

    # get share of residential and sevice floor area
    sec_w = (area_tot / area_tot.groupby(["country"]).sum())["value"]
    # get the total cost-energy-savings weight by sector area
    tot = dE_and_cost.apply(lambda col: col * sec_w, axis=0).groupby(level=0).sum()
    tot.set_index(pd.MultiIndex.from_product([list(tot.index), ["tot"]]),
                  inplace=True)
    dE_and_cost = dE_and_cost.append(tot).unstack().stack()

    summed_area = pd.DataFrame(area_tot.groupby("country").sum())
    summed_area.set_index(pd.MultiIndex.from_product(
                          [list(summed_area.index), ["tot"]]), inplace=True)
    area_tot = area_tot.append(summed_area).unstack().stack()

# %% ******* (4) SAVE   ************************************************

    dE_and_cost.to_csv(snakemake.output.retro_cost)
    area_tot.to_csv(snakemake.output.floor_area)



