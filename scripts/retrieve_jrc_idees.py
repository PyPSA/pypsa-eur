# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract JRC IDEES 2021 data.
"""

import logging
import os
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

# Define the base URL
url_jrc = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/JRC-IDEES/JRC-IDEES-2021_v1/JRC-IDEES-2021.zip"

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_jrc_idees")
        rootpath = ".."
    else:
        rootpath = "."

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    to_fn = snakemake.output[0]
    to_fn_zp = to_fn + ".zip"

    # download .zip file
    logger.info(f"Downloading JRC IDEES from '{url_jrc}'.")
    progress_retrieve(url_jrc, to_fn_zp, disable=disable_progress)

    # extract
    logger.info("Extracting JRC IDEES data.")
    with zipfile.ZipFile(to_fn_zp, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    logger.info(f"JRC IDEES data available in '{to_fn}'.")
