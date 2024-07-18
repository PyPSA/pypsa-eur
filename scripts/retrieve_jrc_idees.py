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

import requests
from _helpers import configure_logging, progress_retrieve, set_scenario_config
from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)

# Define the base URL
url_jrc = (
    "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/JRC-IDEES/JRC-IDEES-2021_v1/"
)

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

    # create a local directory to save the zip files
    local_dir = snakemake.output[0]
    if not os.path.exists(local_dir):
        os.makedirs(local_dir)

    # get the list of zip files from the JRC URL
    response = requests.get(url_jrc)
    soup = BeautifulSoup(response.text, "html.parser")
    zip_files = [
        link.get("href")
        for link in soup.find_all("a")
        if link.get("href").endswith(".zip")
    ]

    logger.info(
        f"Downloading {len(zip_files)} .zip files for JRC IDEES from '{url_jrc}'."
    )

    # download and unpack each zip file
    for zip_file in zip_files:
        logger.info(f"Downloading and unpacking {zip_file}")
        zip_url = url_jrc + zip_file
        to_fn = local_dir + "/" + zip_file[:-4]
        progress_retrieve(zip_url, to_fn, disable=disable_progress)
