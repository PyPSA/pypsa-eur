# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020 @JanFrederickUnnewehr, The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script calculates the cement production in Europe, assuming that most of it is satisfied locally, thus excluding trade (aroung 6%)
Then the nodal produciton is based on the existing plant capacities.

"""

import logging
import country_converter as coco

logger = logging.getLogger(__name__)
import re
import numpy as np
import pandas as pd
from _helpers import configure_logging

cc = coco.CountryConverter()

# -------------------- CEMENT PRODUCTION ------------------------

def cement_projections(coeffs, investment_years):
    """
    Function to project the cement production based on the planning horizon parameter and node

    """
    ssp_data = pd.read_excel(coeffs)
    gdp = ssp_data[ssp_data["Variable"] == "GDP|PPP"]
    population = ssp_data[ssp_data["Variable"] == "Population"]

    gdp.set_index("Region", inplace=True)
    population.set_index("Region", inplace=True)

    columns2remove = ["Model", "Scenario", "Variable", "Unit", "2020"]
    gdp = gdp.drop(columns=columns2remove, errors="ignore")
    population = population.drop(columns=columns2remove, errors="ignore")

    countries = snakemake.params.countries

    population.index = cc.convert(population.index, to="iso2")
    population = population[population.index.isin(countries)]
    population_eu = population.sum()

    gdp.index = cc.convert(gdp.index, to="iso2")
    gdp = gdp[gdp.index.isin(countries)]
    gdp_eu = gdp.sum()

    gdppc_eu = gdp_eu * 1e3 / population_eu # $ / person

    # Regression parameters for Non Linear of paper van Ruujien et al. 2016
    a = 487
    b = -3047
    e_component = a * np.exp(b / gdppc_eu)
    cement_prod_eu = (e_component * population_eu) # kt/yr

    # Values are not too in line with the actual production so we take only the trend and we calculate the change with respect to 2025 for each year
    cement_prod_eu_change = ((cement_prod_eu - cement_prod_eu.iloc[0] ) / cement_prod_eu.iloc[0] )

    # Retrieve historical data from JRC - IDEES and Excel sheet with manually retrieved data
    idees = pd.read_excel(f"{snakemake.input.idees}/EU27/JRC-IDEES-2021_Industry_EU27.xlsx",sheet_name='NMM',index_col=0,header=0,)
    idees_2021 = idees.loc["Cement (kt)", 2021].iloc[0]
    cement_prod_extra_eu = pd.read_excel(f"{snakemake.input.cement_extra_eu}", index_col=0,header=0)
    cement_prod_extra_eu = cement_prod_extra_eu.sum()['kt/yr']
    tot_cement_prod = idees_2021 + cement_prod_extra_eu
    future_eu_cement_production = tot_cement_prod * (1 + cement_prod_eu_change)

    # Retrieve cement plants to distribute the future european cement production per country
    cement_plants = pd.read_excel(f"{snakemake.input.cement_plants}", sheet_name="SFI_ALD_Cement_Database", index_col=0, header=0)
    cement_plants.loc[:,'country'] = cc.convert(cement_plants.loc[:,'country'], to="ISO2")
    cement_plants = cement_plants[cement_plants['country'].isin(countries)]
    keys = pd.read_csv(snakemake.input.industrial_distribution_key, index_col=0)
    keys = keys.loc[:,'Cement']
    country_capacity = (cement_plants[['country', 'capacity']].groupby('country').sum())
    node_capacity = pd.DataFrame(0,index=keys.index, columns=['capacity'])

    index = node_capacity.index
    node_capacity.index = node_capacity.index.str[:2] 
    country_capacity.index = country_capacity.index.str[:2]
    node_capacity['capacity'] = node_capacity.index.map(country_capacity['capacity'])
    node_capacity.index = index

    node_capacity['capacity'] = node_capacity['capacity'] * keys
    node_capacity['capacity'] = node_capacity['capacity'].fillna(0)
    total_capacity = node_capacity['capacity'].sum()  
    node_capacity['share'] = node_capacity['capacity'] / total_capacity

    #nodal_cement_prod = pd.DataFrame(0,index=node_capacity.index, columns=future_eu_cement_production.index)
    nodal_cement_prod = pd.DataFrame(
        np.outer(node_capacity['share'], future_eu_cement_production),
        index=node_capacity.index,
        columns=future_eu_cement_production.index
    )
    
    return nodal_cement_prod



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_steel_production_projections",
            simpl="",
            clusters="39",
            planning_horizons=2050,
            run="baseline",
        )

    configure_logging(snakemake)

    investment_years = [2025, 2030, 2035, 2040, 2045, 2050]
    ssp = snakemake.input.ssp

    nodal_cement_prod = cement_projections(ssp, investment_years)
    nodal_cement_prod.to_csv(snakemake.output.cement_prod)
