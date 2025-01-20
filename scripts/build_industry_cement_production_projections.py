# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020 @JanFrederickUnnewehr, The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
This script calculates cement production in Europe under the assumption that most production is satisfied locally,
with minimal trade (~6%). Future cement production is projected based on GDP and population trends, 
using regression parameters from the study by van Ruijven et al. (2016). 
The nodal distribution of cement production is determined based on existing plant capacities.

The script provides three policy scenarios for European cement production:
1. **Maintain Quota**: Europe continues producing cement at a fixed percentage of global demand, ensuring stable local production.
2. **Reindustrialization**: Europe increases its cement production, reflecting policies that promote industrial growth and self-reliance.
3. **Deindustrialization**: Europe’s cement production follows a downward trend, reflecting continued industrial decline and reduced local production.
"""

import logging
import country_converter as coco

logger = logging.getLogger(__name__)
import numpy as np
import pandas as pd
from _helpers import configure_logging

cc = coco.CountryConverter()

# -------------------- CEMENT PRODUCTION ------------------------

def cement_projections(coeffs, investment_years):
    """
    Projects European cement production based on GDP and population trends and distributes production across nodes
    using existing plant capacities.

    Parameters:
    - coeffs: str, path to the input file with GDP and population data.
    - investment_years: list of int, years for which cement production is projected.

    Returns:
    - nodal_cement_prod: pd.DataFrame, nodal cement production over the projection horizon.
    """
    # Load GDP and population data
    ssp_data = pd.read_excel(coeffs)
    gdp = ssp_data[ssp_data["Variable"] == "GDP|PPP"]
    population = ssp_data[ssp_data["Variable"] == "Population"]

    # Set country identifiers as index for GDP and population data
    gdp.set_index("Region", inplace=True)
    population.set_index("Region", inplace=True)

    # Remove unnecessary columns
    columns2remove = ["Model", "Scenario", "Variable", "Unit", "2020"]
    gdp = gdp.drop(columns=columns2remove, errors="ignore")
    population = population.drop(columns=columns2remove, errors="ignore")

    # Filter data for modeled European countries
    countries = snakemake.params.countries
    population.index = cc.convert(population.index, to="iso2")
    population = population[population.index.isin(countries)]
    population_eu = population.sum()

    gdp.index = cc.convert(gdp.index, to="iso2")
    gdp = gdp[gdp.index.isin(countries)]
    gdp_eu = gdp.sum()

    # Calculate GDP per capita (2017 USD per person)
    gdppc_eu = gdp_eu * 1e3 / population_eu

    # Use regression parameters from van Ruijven et al. (2016) to estimate cement production
    a = 487  # Scaling factor
    b = -3047  # Exponential parameter
    e_component = a * np.exp(b / gdppc_eu)
    cement_prod_eu = e_component * population_eu  # kt/year

    # Adjust production trends using 2025 as a baseline
    cement_prod_eu_change = (cement_prod_eu - cement_prod_eu.iloc[0]) / cement_prod_eu.iloc[0]

    # Load historical production data (2021 as the base year)
    idees = pd.read_excel(
        f"{snakemake.input.idees}/EU27/JRC-IDEES-2021_Industry_EU27.xlsx",
        sheet_name='NMM',
        index_col=0,
        header=0,
    )
    idees_2021 = idees.loc["Cement (kt)", 2021].iloc[0]

    # Include extra-European cement production
    cement_prod_extra_eu = pd.read_excel(f"{snakemake.input.cement_extra_eu}", index_col=0, header=0)
    cement_prod_extra_eu = cement_prod_extra_eu.sum()['kt/yr']
    tot_cement_prod = idees_2021 + cement_prod_extra_eu

    # Project future European cement production
    future_eu_cement_production = tot_cement_prod * (1 + cement_prod_eu_change)

    # Load existing cement plant data to distribute production
    cement_plants = pd.read_excel(
        f"{snakemake.input.cement_plants}",
        sheet_name="SFI_ALD_Cement_Database",
        index_col=0,
        header=0,
    )
    cement_plants.loc[:, 'country'] = cc.convert(cement_plants.loc[:, 'country'], to="ISO2")
    cement_plants = cement_plants[cement_plants['country'].isin(countries)]

    # Load industrial distribution keys and aggregate country capacities
    keys = pd.read_csv(snakemake.input.industrial_distribution_key, index_col=0)
    keys = keys.loc[:, 'Cement']
    country_capacity = cement_plants[['country', 'capacity']].groupby('country').sum()
    node_capacity = pd.DataFrame(0, index=keys.index, columns=['capacity'])

    # Map country capacities to nodes
    index = node_capacity.index
    node_capacity.index = node_capacity.index.str[:2]
    country_capacity.index = country_capacity.index.str[:2]
    node_capacity['capacity'] = node_capacity.index.map(country_capacity['capacity'])
    node_capacity.index = index

    # Adjust nodal capacity by distribution keys
    node_capacity['capacity'] = node_capacity['capacity'] * keys
    node_capacity['capacity'] = node_capacity['capacity'].fillna(0)

    # Calculate nodal share of total capacity
    total_capacity = node_capacity['capacity'].sum()
    node_capacity['share'] = node_capacity['capacity'] / total_capacity

    # Distribute projected cement production to nodes
    nodal_cement_prod = pd.DataFrame(
        np.outer(node_capacity['share'], future_eu_cement_production),
        index=node_capacity.index,
        columns=future_eu_cement_production.index,
    )

    return nodal_cement_prod


if __name__ == "__main__":
    # Mock snakemake if running script standalone
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_steel_production_projections",
            simpl="",
            clusters="39",
            planning_horizons=2050,
            run="baseline",
        )

    # Configure logging
    configure_logging(snakemake)

    # Define investment years for projections
    investment_years = [2025, 2030, 2035, 2040, 2045, 2050]
    ssp = snakemake.input.ssp

    # Project and distribute nodal cement production
    nodal_cement_prod = cement_projections(ssp, investment_years)
    nodal_cement_prod.to_csv(snakemake.output.cement_prod)
