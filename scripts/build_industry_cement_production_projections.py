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


# Define the function to calculate the cement production per capita [kg/cap/yr] based on https://www.sciencedirect.com/science/article/pii/S0921344916301008?via%3Dihub#bibl0005 
def calculate_eu_cement_production(gdppc, population, investment_years):

    # Create the output dataframe
    cement_demand = pd.DataFrame(index=gdppc.index, columns=[investment_years])

    # Regression parameters for Non Linear
    a = 487
    b = -3047

    # Calculate the e^x component for all countries and years
    e_component = a * np.exp(b / gdppc)
    e_component.columns = e_component.columns.astype(str)
    cement_demand = (e_component * population) # kt/yr
    
    cement_demand.index = cc.convert(cement_demand.index, to="iso2")
    #cement_demand = cement_demand[cement_demand.index.isin(countries)]
    print(f"European EU cement demand {cement_demand.sum()}")
     # I need to sum first the GDP and then divide by the population for europe

    return cement_demand


def cement_projections(coeffs, investment_years):
    """
    Function to project the cement production based on the planning horizon parameter and country

    """

    ssp_data = pd.read_excel(coeffs)
    gdp = ssp_data[ssp_data["Variable"] == "GDP|PPP"]
    population = ssp_data[ssp_data["Variable"] == "Population"]

    columns2remove = ["Model", "Scenario", "Variable", "Unit", "2020"]

    # Set 'Region' column as the index
    gdp.set_index("Region", inplace=True)
    population.set_index("Region", inplace=True)

    # Remove specified columns
    gdp = gdp.drop(columns=columns2remove, errors="ignore")
    population = population.drop(columns=columns2remove, errors="ignore")

    # Retrieve the countries selected in the config file
    countries = snakemake.params.countries
    
    population.index = cc.convert(population.index, to="iso2")
    population = population[population.index.isin(countries)]
    population_eu = population.sum()

    gdp.index = cc.convert(gdp.index, to="iso2")
    gdp = gdp[gdp.index.isin(countries)]
    gdp_eu = gdp.sum()

    gdppc_eu = gdp_eu * 1e3 / population_eu # $ / person
        # Create the output dataframe
    #cement_demand = pd.DataFrame(index=gdppc.index, columns=[investment_years])

    # Regression parameters for Non Linear
    a = 487
    b = -3047

    # Calculate the e^x component for all countries and years
    e_component = a * np.exp(b / gdppc_eu)
    #e_component.columns = e_component.columns.astype(str)
    cement_demand_eu = (e_component * population_eu) # kt/yr
    
    print(f"European EU cement demand {cement_demand_eu}")


    # Calculate the percentage change with respect to 2025 for each year
    cement_demand_eu_change = ((cement_demand_eu - cement_demand_eu.iloc[0] ) / cement_demand_eu.iloc[0] ) * 100

    # Retrieve historical data from JRC - IDEES 
    idees = pd.read_excel(f"{snakemake.input.idees}/EU27/JRC-IDEES-2021_Industry_EU27.xlsx",sheet_name='NMM',index_col=0,header=0,)
    idees_2021 = idees.loc["Cement (kt)", 2021][0]
    albania_2024 = 2.8 * 1e3 # kt/yr Source: https://www.globalcement.com/news/item/17798-update-on-the-central-balkans-august-2024 
    bosnia_2024 = 1.6 * 1e3 #kt/yr Source: https://www.globalcement.com/news/item/17798-update-on-the-central-balkans-august-2024 
    montenegro_2024 = 0 #kt/yr Source: https://www.globalcement.com/news/item/17798-update-on-the-central-balkans-august-2024 
    normace_2024 = 1.4 * 1e3 #kt/yr Source: https://www.globalcement.com/news/item/17798-update-on-the-central-balkans-august-2024 
    serbia_2024 = 2.7 * 1e3 #kt/yr Source: https://www.globalcement.com/news/item/17798-update-on-the-central-balkans-august-2024 
    kosovo_2024 = 0.5 * 1e3 #kt/yr Source: https://www.globalcement.com/news/item/17798-update-on-the-central-balkans-august-2024 
    uk_2021 = 7.4 * 1e3 #kt/yr Source: https://gmk.center/en/news/uk-reduced-steel-production-by-6-5-y-y-in-2023/#:~:text=Steel%20production%20in%20the%20country,24th%20among%20steel%20producing%20countries. 
    norway_2021 = 1.7 * 1e3 # kt/yr Source: https://www.globalcement.com/magazine/articles/1199-cement-in-northern-europe + https://www.sement.heidelbergmaterials.no/en/NorcemKjopsvik_en 
    switzerland_2021 = 4.87 * 1e3 # kt/yr Source: https://www.sciencedirect.com/science/article/pii/S0959652620354597
    
    return cement_prod


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

    # Selected years for the projections
    investment_years = [2025, 2030, 2035, 2040, 2045, 2050]

    ssp = snakemake.input.ssp

    # Project the load based on empirical analyses for future years scenarios

    eu_prod = cement_projections(ssp, investment_years)
    eu_prod.to_csv(snakemake.output.cement_demand)
