# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2021-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve gas infrastructure data from
https://zenodo.org/record/4767098/files/IGGIELGN.zip.
"""

import logging
import zipfile
from pathlib import Path

from _helpers import progress_retrieve

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_gas_network_data")
        rootpath = ".."
    else:
        rootpath = "."

    url = "https://zenodo.org/record/4767098/files/IGGIELGN.zip"

    # Save locations
    zip_fn = Path(f"{rootpath}/IGGIELGN.zip")
    to_fn = Path(f"{rootpath}/data/gas_network/scigrid-gas")

    logger.info(f"Downloading databundle from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, zip_fn, disable=disable_progress)

    logger.info("Extracting databundle.")
    zipfile.ZipFile(zip_fn).extractall(to_fn)

    zip_fn.unlink()

    logger.info(f"Gas infrastructure data available in '{to_fn}'.")
