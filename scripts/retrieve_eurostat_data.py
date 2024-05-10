# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract eurostat energy balances data.
"""


import logging
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_eurostat_data")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    url_eurostat = "https://ec.europa.eu/eurostat/documents/38154/4956218/Balances-December2022.zip/f7cf0d19-5c0f-60ad-4e48-098a5ddd6e48?t=1671184070589"
    tarball_fn = Path(f"{rootpath}/data/eurostat/eurostat_2023.zip")
    to_fn = Path(
        f"{rootpath}/data/eurostat/eurostat-energy_balances-april_2023_edition/"
    )

    logger.info(f"Downloading Eurostat data from '{url_eurostat}'.")
    progress_retrieve(url_eurostat, tarball_fn, disable=disable_progress)

    logger.info("Extracting Eurostat data.")
    with zipfile.ZipFile(tarball_fn, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    logger.info(f"Eurostat data available in '{to_fn}'.")
