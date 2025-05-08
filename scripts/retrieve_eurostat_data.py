# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
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
    url_eurostat = (
        # "https://ec.europa.eu/eurostat/documents/38154/4956218/Balances-April2023.zip" # link down
        "https://tubcloud.tu-berlin.de/s/prkJpL7B9M3cDPb/download/Balances-April2023.zip"
    )
    tarball_fn = Path(f"{rootpath}/data/eurostat/eurostat_2023.zip")
    to_fn = Path(f"{rootpath}/data/eurostat/Balances-April2023/")

    logger.info(f"Downloading Eurostat data from '{url_eurostat}'.")
    progress_retrieve(url_eurostat, tarball_fn, disable=disable_progress)

    logger.info("Extracting Eurostat data.")
    with zipfile.ZipFile(tarball_fn, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    logger.info(f"Eurostat data available in '{to_fn}'.")
