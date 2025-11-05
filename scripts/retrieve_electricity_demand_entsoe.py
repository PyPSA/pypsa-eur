# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve electricity demand data from ENTSOE.
"""

import logging

import pandas as pd
from entsoe import EntsoePandasClient

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

# avoid unnecessary API calls by once retrieving all supported countries
COUNTRIES = [
    "AL",
    "AT",
    "BE",
    "BA",
    "BG",
    "CH",
    "CY",
    "CZ",
    "DE",
    "DK",
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

        snakemake = mock_snakemake("retrieve_electricity_demand_entsoe")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    token = snakemake.params.entsoe_token

    if token:
        logger.info(
            "ENTSOE token provided. Retrieving live data from transparency platform."
        )
        client = EntsoePandasClient(api_key=token)
        start = pd.Timestamp("20150101", tz="UTC")
        end = pd.Timestamp.today(tz="UTC")
        loads = []
        for country_code in COUNTRIES:
            logger.info(f"Querying load for {country_code}...")
            loads.append(
                client.query_load(country_code, start=start, end=end)
                .squeeze()
                .rename(country_code)
                .tz_convert("UTC")
                .resample("1h")
                .mean()
            )

        df = pd.concat(loads, axis=1, join="outer")

    else:
        logger.info(
            "No ENTSOE token provided. Retrieving pre-built data from data bundle."
        )
        # TODO add to zenodo data bundle once vetted
        prebuilt_url = "https://tubcloud.tu-berlin.de/s/qKKdAiNeDxFscDH/download/electricity_demand_entsoe_raw.csv"
        df = pd.read_csv(prebuilt_url)

    df.to_csv(snakemake.output[0])
