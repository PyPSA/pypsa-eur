# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2021-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract data bundle for sector-coupled studies.
"""

import logging

logger = logging.getLogger(__name__)

import os
import sys
import tarfile
from pathlib import Path

# Add pypsa-eur scripts to path for import of _helpers
sys.path.insert(0, os.getcwd() + "/../pypsa-eur/scripts")

from _helpers import configure_logging, progress_retrieve

if __name__ == "__main__":
    configure_logging(snakemake)

    url = "https://zenodo.org/record/5824485/files/pypsa-eur-sec-data-bundle.tar.gz"

    tarball_fn = Path("sector-bundle.tar.gz")
    to_fn = Path("data")

    logger.info(f"Downloading databundle from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, tarball_fn, disable=disable_progress)

    logger.info("Extracting databundle.")
    tarfile.open(tarball_fn).extractall(to_fn)

    tarball_fn.unlink()

    logger.info(f"Databundle available in '{to_fn}'.")
