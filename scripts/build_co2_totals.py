# SPDX-FileCopyrightText: : 2025 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Calculate historical CO2 emissions per country using EEA and Eurostat data.

Outputs
-------
- ``resources/<run_name>/co2_totals.csv``: CO2 emissions per country and sector.
"""

import logging

import country_converter as coco
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

cc = coco.CountryConverter()
logger = logging.getLogger(__name__)
idx = pd.IndexSlice

eu28 = cc.EU28as("ISO2").ISO2.tolist()
eu28_eea = eu28.copy()
eu28_eea.remove("GB")
eu28_eea.append("UK")

to_ipcc = {
    "electricity": "1.A.1.a - Public Electricity and Heat Production",
    "residential non-elec": "1.A.4.b - Residential",
    "services non-elec": "1.A.4.a - Commercial/Institutional",
    "rail non-elec": "1.A.3.c - Railways",
    "road non-elec": "1.A.3.b - Road Transportation",
    "domestic navigation": "1.A.3.d - Domestic Navigation",
    "international navigation": "1.D.1.b - International Navigation",
    "domestic aviation": "1.A.3.a - Domestic Aviation",
    "international aviation": "1.D.1.a - International Aviation",
    "total energy": "1 - Energy",
    "industrial processes": "2 - Industrial Processes and Product Use",
    "agriculture": "3 - Agriculture",
    "agriculture, forestry and fishing": "1.A.4.c - Agriculture/Forestry/Fishing",
    "LULUCF": "4 - Land Use, Land-Use Change and Forestry",
    "waste management": "5 - Waste management",
    "other": "6 - Other Sector",
    "indirect": "ind_CO2 - Indirect CO2",
    "total wL": "Total (with LULUCF)",
    "total woL": "Total (without LULUCF)",
}


def reverse(dictionary: dict) -> dict:
    return {v: k for k, v in dictionary.items()}


def build_eea_co2(
    input_co2: str, year: int = 1990, emissions_scope: str = "CO2"
) -> pd.DataFrame:
    """
    Calculate CO2 emissions for a given year based on EEA data in Mt.

    Parameters
    ----------
    input_co2 : str
        Path to the input CSV file with CO2 data.
    year : int, optional
        Year for which to calculate emissions, by default 1990.
    emissions_scope : str, optional
        Scope of the emissions to consider, by default "CO2".

    Returns
    -------
    pd.DataFrame
        DataFrame with CO2 emissions for the given year.

    Notes
    -----
    - The function reads the `input_co2` data and for a specific `year` and `emission scope`
    - It calculates "industrial non-elec" and "agriculture" emissions from that data
    - It drops unneeded columns and converts the emissions to Mt.

    References
    ----------
    - `EEA CO2 data <https://www.eea.europa.eu/data-and-maps/data/national-emissions-reported-to-the-unfccc-and-to-the-eu-greenhouse-gas-monitoring-mechanism-16>`_ (downloaded 201228, modified by EEA last on 201221)
    """
    df = pd.read_csv(input_co2, encoding="latin-1", low_memory=False)

    df.replace(dict(Year="1985-1987"), 1986, inplace=True)
    df.Year = df.Year.astype(int)
    index_col = ["Country_code", "Pollutant_name", "Year", "Sector_name"]
    df = df.set_index(index_col).sort_index()

    cts = ["CH", "EUA", "NO"] + eu28_eea

    slicer = idx[cts, emissions_scope, year, to_ipcc.values()]
    emissions = (
        df.loc[slicer, "emissions"]
        .unstack("Sector_name")
        .rename(columns=reverse(to_ipcc))
        .droplevel([1, 2])
    )

    emissions.rename(index={"EUA": "EU28", "UK": "GB"}, inplace=True)

    to_subtract = [
        "electricity",
        "services non-elec",
        "residential non-elec",
        "road non-elec",
        "rail non-elec",
        "domestic aviation",
        "international aviation",
        "domestic navigation",
        "international navigation",
        "agriculture, forestry and fishing",
    ]
    emissions["industrial non-elec"] = emissions["total energy"] - emissions[
        to_subtract
    ].sum(axis=1)

    emissions["agriculture"] += emissions["agriculture, forestry and fishing"]

    to_drop = [
        "total energy",
        "total wL",
        "total woL",
        "agriculture, forestry and fishing",
    ]
    emissions.drop(columns=to_drop, inplace=True)

    # convert from Gt to Mt
    return emissions / 1e3


def build_eurostat_co2(eurostat: pd.DataFrame, year: int = 1990) -> pd.Series:
    """
    Calculate CO2 emissions for a given year based on Eurostat fuel consumption
    data and fuel-specific emissions.

    Parameters
    ----------
    eurostat : pd.DataFrame
        DataFrame with Eurostat data.
    year : int, optional
        Year for which to calculate emissions, by default 1990.

    Returns
    -------
    pd.Series
        Series with CO2 emissions for the given year.

    Notes
    -----
    - The function hard-sets fuel-specific emissions:
        - solid fuels: 0.36 tCO2_equi/MW_th (approximates coal)
        - oil: 0.285 tCO2_equi/MW_th (average of distillate and residue)
        - natural gas: 0.2 tCO2_equi/MW_th
    - It then multiplies the Eurostat fuel consumption data for `year` by the specific emissions and sums the result.

    References
    ----------
    - Oil values from `EIA <https://www.eia.gov/tools/faqs/faq.cfm?id=74&t=11>`_
    - Distillate oil (No. 2)  0.276
    - Residual oil (No. 6)  0.298
    - `EIA Electricity Annual <https://www.eia.gov/electricity/annual/html/epa_a_03.html>`_
    """
    emissions = pd.Series(
        {
            "C0000X0350-0370": 0.36,  # solid fossil fuels
            "O4000XBIO": 0.285,  # oil and petroleum products
            "G3000": 0.2,  # natural gas
        }
    )
    return (
        eurostat.query("year == @year and siec in @emissions.index")
        .assign(value=lambda df: df["value"] * df["siec"].map(emissions))
        .groupby(["country", "nrg_bal"])["value"]
        .sum(min_count=1)
    )


def build_co2_totals(
    countries: list[str], eea_co2: pd.DataFrame, eurostat_co2: pd.DataFrame
) -> pd.DataFrame:
    """
    Combine CO2 emissions data from EEA and Eurostat for a list of countries.

    Parameters
    ----------
    countries : list[str]
        List of country codes for which CO2 totals need to be built.
    eea_co2 : pd.DataFrame
        DataFrame with EEA CO2 emissions data.
    eurostat_co2 : pd.DataFrame
        DataFrame with Eurostat CO2 emissions data.

    Returns
    -------
    pd.DataFrame
        Combined CO2 emissions data for the given countries.

    Notes
    -----
    - The function combines the CO2 emissions from EEA and Eurostat into a single DataFrame for the given countries.
    """
    co2 = eea_co2.reindex(countries)

    for ct in pd.Index(countries).intersection(
        ["BA", "RS", "XK", "AL", "ME", "MK", "UA", "MD"]
    ):
        mappings = {
            "electricity": "TI_EHG_E",
            "residential non-elec": "FC_OTH_HH_E",
            "services non-elec": "FC_OTH_CP_E",
            "road non-elec": "FC_TRA_ROAD_E",
            "rail non-elec": "FC_TRA_RAIL_E",
            "domestic navigation": "FC_TRA_DNAVI_E",
            "international navigation": "INTMARB",
            "domestic aviation": "FC_TRA_DAVI_E",
            "international aviation": "INTAVI",
            # does not include industrial process emissions or fuel processing/refining
            "industrial non-elec": "FC_IND_E",
            # does not include non-energy emissions
            "agriculture": ["FC_OTH_AF_E", "FC_OTH_FISH_E"],
        }

        for i, mi in mappings.items():
            co2.at[ct, i] = eurostat_co2.loc[ct, mi].sum()

    return co2


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_co2_totals")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params.energy
    countries = snakemake.params.countries
    base_year_emissions = params["base_emissions_year"]
    emissions_scope = params["emissions"]

    eurostat = pd.read_csv(snakemake.input.eurostat)
    eea_co2 = build_eea_co2(snakemake.input.co2, base_year_emissions, emissions_scope)
    eurostat_co2 = build_eurostat_co2(eurostat, base_year_emissions)

    co2 = build_co2_totals(countries, eea_co2, eurostat_co2)
    co2.to_csv(snakemake.output.co2_totals)
