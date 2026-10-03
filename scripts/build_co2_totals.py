# SPDX-FileCopyrightText: : 2025 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build historical CO2 or greenhouse gas emissions per country, year and sector from the EDGAR inventory.

The EDGAR emissions by country and IPCC 2006 source category are mapped to the
sectors of PyPSA-Eur. Emissions of the energy sector not mapped to a specific
sector (e.g. refineries, manufacturing, fugitive emissions) are assigned to
"industrial non-elec". EDGAR reports residential, services and agricultural
fuel combustion (1.A.4) as one sector and international aviation and shipping
only as global totals. These bunkers are therefore not included. EDGAR reports
Serbia, Montenegro and Kosovo as one entity, which is split by population.

Outputs
-------
- ``resources/<run_name>/co2_totals.csv``: Emissions in Mt per country, year and sector.

References
----------
- `EDGAR Community GHG Database <https://edgar.jrc.ec.europa.eu/dataset_ghg2026>`_
- Crippa, M. et al. (2026), GHG emissions of all world countries - 2026 Report, doi:10.2760/7717504
"""

import logging

import country_converter as coco
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

# Longest matching prefix of the IPCC 2006 code determines the sector
IPCC_SECTORS = {
    "1": "industrial non-elec",
    "1.A.1.a": "electricity",
    "1.A.3.a": "domestic aviation",
    "1.A.3.b": "road non-elec",
    "1.A.3.c": "rail non-elec",
    "1.A.3.d": "domestic navigation",
    "1.A.4": "buildings non-elec",
    "2": "industrial processes",
    "3": "agriculture",
    "4": "waste management",
    "5.A": "indirect",
    "5.B": "industrial non-elec",  # fossil fuel fires
}

# Population in 2020 from World Bank (SP.POP.TOTL) to split "Serbia and Montenegro"
SCG_POPULATION = {"RS": 6899126, "ME": 626590, "XK": 1790151}


def ipcc_to_sector(code: str) -> str:
    """
    Map an IPCC 2006 source category code to a PyPSA-Eur emission sector.

    Parameters
    ----------
    code : str
        IPCC 2006 code, e.g. "1.A.3.b".

    Returns
    -------
    str
        Emission sector.
    """
    matches = [k for k in IPCC_SECTORS if code == k or code.startswith(k + ".")]
    if not matches:
        raise ValueError(f"IPCC code '{code}' is not mapped to an emission sector.")
    return IPCC_SECTORS[max(matches, key=len)]


def build_co2_totals(fn: str, countries: list[str]) -> pd.DataFrame:
    """
    Read EDGAR emissions and aggregate them by country, year and sector.

    Parameters
    ----------
    fn : str
        Path to the EDGAR Excel file with emissions by country and IPCC 2006 code.
    countries : list[str]
        ISO2 country codes to include.

    Returns
    -------
    pd.DataFrame
        Emissions in Mt with a (country, year) index and sectors as columns.
    """
    df = pd.read_excel(fn, sheet_name="IPCC 2006", skiprows=9)
    years = df.filter(like="Y_").columns

    codes = df.Country_code_A3.unique()
    iso2 = dict(zip(codes, coco.convert(codes, to="ISO2", not_found=None)))
    df["country"] = df.Country_code_A3.map(iso2)

    scg = df.query("Country_code_A3 == 'SCG'")
    shares = pd.Series(SCG_POPULATION) / sum(SCG_POPULATION.values())
    df = pd.concat(
        [df.query("Country_code_A3 != 'SCG'")]
        + [scg.assign(country=ct, **scg[years].mul(s)) for ct, s in shares.items()]
    )

    missing = pd.Index(countries).difference(df.country.unique())
    if not missing.empty:
        raise ValueError(f"No EDGAR emissions for countries {missing.tolist()}.")

    emissions = (
        df.query("country in @countries")
        .assign(
            sector=lambda d: d.ipcc_code_2006_for_standard_report.map(ipcc_to_sector)
        )
        .melt(id_vars=["country", "sector"], value_vars=years, var_name="year")
        .assign(year=lambda d: d.year.str.removeprefix("Y_").astype(int))
        .pivot_table(
            index=["country", "year"],
            columns="sector",
            values="value",
            aggfunc="sum",
            fill_value=0.0,
        )
    )

    # convert from Gg to Mt
    return emissions / 1e3


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_co2_totals")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    co2 = build_co2_totals(snakemake.input.edgar, snakemake.params.countries)
    co2.to_csv(snakemake.output.co2_totals)
