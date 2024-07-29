# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve ammonia demand from https://www.usgs.gov/centers/national-minerals-information-center/nitrogen-statistics-and-information.
"""

import logging
import os
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

# Define the base URL
url = "https://d9-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/s3fs-public/media/files/myb1-2022-nitro-ert.xlsx"

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_ammonia_demand")
        rootpath = ".."
    else:
        rootpath = "."

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    to_fn = snakemake.output[0]


    # download .zip file
    logger.info(f"Downloading Ammonia demand from {url}.")
    progress_retrieve(url, to_fn, disable=disable_progress)


    logger.info(f"Ammonia demand data available in '{to_fn}'.")
