# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2021-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract data bundle for sector-coupled studies.
"""

import logging

logger = logging.getLogger(__name__)

import tarfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_databundle")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)

    url = "https://zenodo.org/record/5824485/files/pypsa-eur-sec-data-bundle.tar.gz"

    tarball_fn = Path(f"{rootpath}/sector-bundle.tar.gz")
    to_fn = Path(rootpath) / Path(snakemake.output[0]).parent.parent

    logger.info(f"Downloading databundle from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, tarball_fn, disable=disable_progress)

    logger.info("Extracting databundle.")
    tarfile.open(tarball_fn).extractall(to_fn)

    tarball_fn.unlink()

    logger.info(f"Databundle available in '{to_fn}'.")
