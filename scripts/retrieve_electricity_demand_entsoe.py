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

    assert token is not None, "ENTSOE API token must be provided!"

    client = EntsoePandasClient(api_key=token)
    start = pd.Timestamp("20150101", tz="UTC")
    end = pd.Timestamp.today(tz="UTC")
    country_code = snakemake.wildcards.country

    logger.info(f"Querying load for {country_code}...")

    df = (
        client.query_load(country_code, start=start, end=end)
        .squeeze()
        .rename(country_code)
        .tz_convert("UTC")
        .resample("1h")
        .mean()
    )

    df.to_csv(snakemake.output.csv)
