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

    # Retrieve the countries selected in the config file
    countries = snakemake.params.countries

    for year in investment_years:

        a = 487
        b = -3047

        cement_demand.loc[country, year] = (
            (
                a*np.e**(b/gdppc[year][country])
            )
            * population[year][country] #milion but we then divide by 1e6 to get kt from kg
        )  # kt/yr

    cement_demand.index = cc.convert(cement_demand.index, to="iso2")
    cement_demand = cement_demand[cement_demand.index.isin(countries)]

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
    gdppc = gdp * 1000 / population  # 2017USD/person

    # Convert string column names containing numbers to numeric values
    gdppc.columns = pd.to_numeric(gdppc.columns, errors="coerce")
    gdp.columns = pd.to_numeric(gdp.columns, errors="coerce")

    cement_prod = calculate_eu_cement_production(
        gdppc, population, investment_years
    )  # ktons of cement 

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
