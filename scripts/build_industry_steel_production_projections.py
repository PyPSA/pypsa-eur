# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020 @JanFrederickUnnewehr, The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
This script calculates global steel demand using the intensity of use (IOU) hypothesis developed by Linda Wårell.
The hypothesis identifies GDP per capita as the primary driver of steel consumption, establishing a relationship between economic development and material demand.

The script estimates global steel demand for selected investment years and refines this analysis to focus on Europe’s contribution to global production, specifically for countries modeled within PyPSA-Eur.
After determining global steel demand, the script provides three possible scenarios for European production:

- Maintain Quota: Europe continues to produce a fixed percentage of global steel demand, maintaining constant levels of internal consumption and export.
- Reindustrialization: Europe reverses its deindustrialization trend by producing at the maximum historical share of global steel demand, reflecting aggressive reindustrialization policies.
- Deindustrialization: Europe follows its current trajectory of deindustrialization, with its share of global production decreasing linearly over time.

This approach enables a flexible examination of Europe’s steel production under varying policy and economic conditions.

"""

import logging
import re
import dateutil
import numpy as np
import pandas as pd
from _helpers import configure_logging
from pandas import Timedelta as Delta

# Setup logging
logger = logging.getLogger(__name__)

# -------------------- STEEL PRODUCTION ------------------------

# Define the function to calculate the intensity of use for steel
def calculate_eu_steel_production(gdppc, gdp, investment_years):
    """
    Calculate the projected steel demand in Europe using the IOU approach.
    
    The IOU method estimates material demand based on GDP per capita (gdppc)
    over time, considering nonlinear relationships with specific coefficients 
    (alpha0, alpha1, alpha2, and alpha3) derived from historical data, as discussed in the referenced paper.

    Parameters:
    - gdppc: DataFrame containing GDP per capita values for various countries and years.
    - investment_years: List of years for which projections are calculated.

    Returns:
    - eu_prod: Series containing projected steel production in Europe (in kilotons, kt).
    """

    # Create an output DataFrame to store projected steel demand per country and year
    steel_demand = pd.DataFrame(index=gdppc.index, columns=[investment_years])

    # Loop over investment years to compute steel demand
    for year in investment_years:
        # Time variable t (1 in 1970, 56 in 2025, etc.)
        t = year - 1970 + 1

        # Loop over each country to calculate steel demand
        for country in gdppc.index:

            # Determine coefficients (alpha0, alpha1, alpha2, alpha3) based on GDP per capita thresholds
            # The thresholds and coefficients are informed by the paper's empirical results
            if gdppc.loc[country, 2025] <= 6000:  # Low-income countries
                alpha0 = 0.044
                alpha1 = -1.5 * 1e-7
                alpha2 = -1.48 * 1e-10
                alpha3 = -0.0004

            elif gdppc.loc[country, 2025] > 6000 and gdppc.loc[country, 2025] <= 11000:  # Middle-income countries
                alpha0 = 0.022
                alpha1 = 2.05 * 1e-6
                alpha2 = -4.22 * 1e-11
                alpha3 = -0.0002

            else:  # High-income countries
                alpha0 = 0.038
                alpha1 = -9.95 * 1e-7
                alpha2 = 8.92 * 1e-12
                alpha3 = -0.0001

            # Calculate steel demand using the IOU formula from the paper
            steel_demand.loc[country, year] = (
                (
                    alpha0
                    + alpha1 * gdppc[year][country]
                    + alpha2 * gdppc[year][country] ** 2
                    + alpha3 * t
                )
                * gdp[year][country]
                * 1000  # Convert from Mt to kt
            )

    # Sum steel demand globally for future world demand
    future_world_demand = steel_demand.sum()

    # Assume Europe maintains a constant share of global production
    # Average share of European production (2019-2023): 8.77% (source: World Steel Association)
    eu_average_prod = 0.0877

    # Calculate Europe's steel production as a share of global demand
    eu_prod = future_world_demand * eu_average_prod

    return eu_prod


def steel_projections(coeffs, investment_years):
    """
    Load SSP (Shared Socioeconomic Pathways) data to project steel demand.
    Data are taken from IIASA SSP Scenario Explorer https://data.ece.iiasa.ac.at/ssp/#/login

    Parameters:
    - coeffs: Path to the SSP data file containing GDP and population projections.
    - investment_years: List of years for which projections are calculated.

    Returns:
    - steel_prod: Projected steel production in Europe (in kilotons, kt).
    """

    # Load SSP data from an Excel file
    ssp_data = pd.read_excel(coeffs)

    # Filter GDP and population data
    gdp = ssp_data[ssp_data["Variable"] == "GDP|PPP"]
    population = ssp_data[ssp_data["Variable"] == "Population"]

    # Columns to exclude from analysis
    columns2remove = ["Model", "Scenario", "Variable", "Unit", "2020"]

    # Set 'Region' as the index for both GDP and population
    gdp.set_index("Region", inplace=True)
    population.set_index("Region", inplace=True)

    # Drop unnecessary columns
    gdp = gdp.drop(columns=columns2remove)
    population = population.drop(columns=columns2remove)

    # Calculate GDP per capita (gdppc)
    # GDP is multiplied by 1000 to adjust units (consistent with SSP datasets)
    gdppc = gdp * 1000 / population  # in 2017 USD/person

    # Convert column names (years) from strings to numeric values
    gdppc.columns = pd.to_numeric(gdppc.columns)
    gdp.columns = pd.to_numeric(gdp.columns)

    # Calculate European steel production using the IOU approach
    steel_prod = calculate_eu_steel_production(
        gdppc, gdp, investment_years
    )  # Output in kt

    return steel_prod


if __name__ == "__main__":
    # If run directly (outside the Snakemake pipeline), use a mock configuration
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        # Mock Snakemake inputs for standalone testing
        snakemake = mock_snakemake(
            "build_industry_steel_production_projections",
            simpl="",
            clusters="39",
            planning_horizons=2050,
            run="baseline",
        )

    # Configure logging for the script
    configure_logging(snakemake)

    # Define the investment years for projections
    investment_years = [2025, 2030, 2035, 2040, 2045, 2050]

    # Path to SSP data input file
    ssp = snakemake.input.ssp

    # Project steel demand and save output
    eu_prod = steel_projections(ssp, investment_years)
    eu_prod.to_csv(snakemake.output.steel_demand)
