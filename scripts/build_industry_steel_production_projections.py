# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020 @JanFrederickUnnewehr, The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script calculates the global steel demand. Then we assume that Europe will continue to produce the same quantity, leaving constant internat consumption and export.

"""

import logging

logger = logging.getLogger(__name__)
import dateutil
import numpy as np
import pandas as pd
import re
from _helpers import configure_logging
from pandas import Timedelta as Delta

# -------------------- STEEL PRODUCTION ------------------------



# Define the function to calculate the intensity of use for steel based on https://www.sciencedirect.com/science/article/pii/S0301420713001207 
def calculate_eu_steel_production(gdppc,gdp,investment_years):

    # Create the output dataframe
    steel_demand = pd.DataFrame(index = gdppc.index, columns = [investment_years])

    for year in investment_years:

        # The time variable t has a value of 1 in 1970, 2 in 1971, so a value of 56 in 2025 and so on
        t = year - 1970 + 1
                
        for country in gdppc.index:

            if gdppc.loc[country, 2025] <= 6000:
                alpha0 = 0.044
                alpha1 = -1.5*1e-7
                alpha2 = -1.48*1e-10
                alpha3 = -0.0004

            elif gdppc.loc[country, 2025] > 6000 and gdppc.loc[country, 2025] <= 11000 :
                alpha0 = 0.022
                alpha1 = 2.05*1e-6
                alpha2 = -4.22*1e-11
                alpha3 = -0.0002
            
            else: 
                alpha0 = 0.038
                alpha1 = -9.95*1e-7
                alpha2 = 8.92*1e-12
                alpha3 = -0.0001

            steel_demand.loc[country, year] = (alpha0 + alpha1 * gdppc[year][country] + alpha2 * gdppc[year][country]**2 + alpha3 * t) * gdp[year][country] *1000 #from Mt to kt
        
    future_world_demand = steel_demand.sum()
    eu_average_prod = 0.0877 # Share of global steel production carried out in Europe, average over 5 years (2019-2023)
    # https://worldsteel.org/steel-topics/statistics/annual-production-steel-data/?ind=P1_crude_steel_total_pub/CHN/IND

    eu_prod = future_world_demand * eu_average_prod

    return eu_prod


def steel_projections(coeffs,investment_years):
    """
    Function to project the steel demand based on the planning horizon parameter and country
    
    """

    ssp_data = pd.read_excel(coeffs)
    gdp = ssp_data[ssp_data['Variable'] == 'GDP|PPP']
    population = ssp_data[ssp_data['Variable'] == 'Population']

    columns2remove = ['Model','Scenario','Variable','Unit','2020']

    # Set 'Region' column as the index
    gdp.set_index('Region', inplace=True)
    population.set_index('Region', inplace=True)

    # Remove specified columns
    gdp = gdp.drop(columns=columns2remove, errors='ignore')
    population = population.drop(columns=columns2remove, errors='ignore')
    gdppc = gdp*1000/population # 2017USD/person

    # Convert string column names containing numbers to numeric values
    gdppc.columns = pd.to_numeric(gdppc.columns, errors='coerce')
    gdp.columns = pd.to_numeric(gdp.columns, errors='coerce')

    steel_prod = calculate_eu_steel_production(gdppc,gdp,investment_years) #ktons of steel products

    return steel_prod



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_steel_production_projections",
            simpl="",
            clusters="5",
            planning_horizons=2050,
            )

    configure_logging(snakemake)

    # Selected years for the projections
    investment_years = [2025,2030,2035,2040,2045,2050]

    ssp = snakemake.input.ssp

    # Project the load based on empirical analyses for future years scenarios

    eu_prod = steel_projections(ssp,investment_years)
    eu_prod.to_csv(snakemake.output[0])
