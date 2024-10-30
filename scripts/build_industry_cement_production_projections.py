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
    cement_prod_eu = (e_component * population_eu) # kt/yr

    # Calculate the change with respect to 2025 for each year
    cement_prod_eu_change = ((cement_prod_eu - cement_prod_eu.iloc[0] ) / cement_prod_eu.iloc[0] )

    # Retrieve historical data from JRC - IDEES 
    idees = pd.read_excel(f"{snakemake.input.idees}/EU27/JRC-IDEES-2021_Industry_EU27.xlsx",sheet_name='NMM',index_col=0,header=0,)
    idees_2021 = idees.loc["Cement (kt)", 2021][0]
    cement_prod_extra_eu = pd.read_excel(f"{snakemake.input.cement_extra_eu}", index_col=0,header=0)
    cement_prod_extra_eu = cement_prod_extra_eu.sum()['kt/yr']
    tot_cement_prod = idees_2021 + cement_prod_extra_eu
    future_eu_cement_production = tot_cement_prod * (1 + cement_prod_eu_change)
    
    return future_eu_cement_production


def retrieve_cement_plants(eu_prod):

    # Distribute the production per region
    cement_plants = pd.read_excel(f"{snakemake.input.cement_plants}", sheet_name="SFI_ALD_Cement_Database", index_col=0, header=0)


    return region_prod, cement_plants

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
    country_prod, cement_plants = cement_plants(eu_prod)
    eu_prod.to_csv(snakemake.output.cement_demand)
