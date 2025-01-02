# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract eurostat GDP and POP data
"""


import logging

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_eurostat_gdp_pop")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    url_eurostat_pop = (
        "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/demo_r_pjanaggr3?format=TSV"
    )


    logger.info(f"Downloading Eurostat population data from '{url_eurostat_pop}'.")
    progress_retrieve(
        url_eurostat_pop, 
        snakemake.output.eurostat_pop, 
        disable=disable_progress
    )
    logger.info(f"Downloaded Eurostat population data to '{snakemake.output.eurostat_pop}'.")
