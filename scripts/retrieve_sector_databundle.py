# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2021-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract data bundle for sector-coupled studies.
"""

import logging
import tarfile
import zipfile
from pathlib import Path

from _helpers import (
    configure_logging,
    progress_retrieve,
    set_scenario_config,
    validate_checksum,
)

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_databundle")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    url = "https://zenodo.org/record/5824485/files/pypsa-eur-sec-data-bundle.tar.gz"

    tarball_fn = Path(f"{rootpath}/sector-bundle.tar.gz")
    to_fn = Path(rootpath) / Path(snakemake.output[0]).parent.parent

    logger.info(f"Downloading databundle from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, tarball_fn, disable=disable_progress)

    validate_checksum(tarball_fn, url)

    logger.info("Extracting databundle.")
    tarfile.open(tarball_fn).extractall(to_fn)

    tarball_fn.unlink()

    logger.info(f"Databundle available in '{to_fn}'.")

    url_eurostat = "https://ec.europa.eu/eurostat/documents/38154/4956218/Balances-December2022.zip/f7cf0d19-5c0f-60ad-4e48-098a5ddd6e48?t=1671184070589"
    tarball_fn = Path(f"{rootpath}/data/bundle-sector/eurostat_2023.zip")
    to_fn = Path(
        f"{rootpath}/data/bundle-sector/eurostat-energy_balances-april_2023_edition/"
    )

    logger.info(f"Downloading Eurostat data from '{url_eurostat}'.")
    progress_retrieve(url_eurostat, tarball_fn, disable=disable_progress)

    logger.info("Extracting Eurostat data.")
    with zipfile.ZipFile(tarball_fn, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    logger.info(f"Eurostat data available in '{to_fn}'.")
