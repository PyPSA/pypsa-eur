# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve power plant outage data from ENTSOE.
"""

import logging

import pandas as pd
from entsoe import EntsoePandasClient
from entsoe.exceptions import NoMatchingDataError

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

COUNTRIES = [
    "AL",
    "AT",
    "BE",
    "BA",
    "BG",
    "CH",
    "CY",
    "CZ",
    "DE_LU",
    "DE_AT_LU",
    "DK_1",
    "DK_2",
    "EE",
    "ES",
    "FI",
    "FR",
    "GB",
    "GR",
    "HR",
    "HU",
    "IE",
    "IT",
    "LT",
    "LU",
    "LV",
    "MD",
    "ME",
    "MK",
    "NL",
    "NO",
    "PL",
    "PT",
    "RO",
    "RS",
    "SE",
    "SI",
    "SK",
    "UA",
    "XK",
]

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_entsoe_outages")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    token = snakemake.params.entsoe_token
    start_year = snakemake.params.start_year
    end_year = snakemake.params.end_year

    if token:
        logger.info(
            "ENTSOE token provided. Retrieving live data from transparency platform."
        )

        client = EntsoePandasClient(api_key=token)
        units = []
        for year in range(start_year, end_year + 1):
            start = f"{year}-01-01"
            end = f"{year + 1}-01-01"
            start = pd.Timestamp(start, tz="UTC")
            end = pd.Timestamp(end, tz="UTC")
            for c in COUNTRIES:
                logger.info(f"Retrieving data for country {c} year {year}")
                try:
                    ugu = client.query_unavailability_of_generation_units(
                        c, start=start, end=end
                    )
                    units.append(ugu)
                except NoMatchingDataError:
                    print(
                        f"No generation data entries found for country {c} in year {year}."
                    )
                try:
                    upu = client.query_unavailability_of_production_units(
                        c, start=start, end=end
                    )
                    units.append(upu)
                except NoMatchingDataError:
                    print(
                        f"No production data entries found for country {c} in year {year}."
                    )

        df = pd.concat(units) if units else pd.DataFrame()

    else:
        logger.info(
            "No ENTSOE token provided. Retrieving pre-built data from data bundle."
        )
        prebuilt_url = ""
        df = pd.read_csv(prebuilt_url)

    df.to_csv(snakemake.output[0])
