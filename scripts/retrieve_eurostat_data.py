# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract eurostat energy balances data.
"""

import logging
import tempfile
import zipfile
from pathlib import Path

from scripts._helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_eurostat_data")
        rootpath = ".."
    else:
        # set to root path of repo or snakemake module
        rootpath = Path(snakemake.output[0]).parent.parent.parent

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    url_eurostat = (
        # "https://ec.europa.eu/eurostat/documents/38154/4956218/Balances-April2023.zip" # link down
        "https://tubcloud.tu-berlin.de/s/prkJpL7B9M3cDPb/download/Balances-April2023.zip"
    )

    to_fn = Path(f"{rootpath}/data/eurostat/Balances-April2023/")

    logger.info(f"Downloading Eurostat data from '{url_eurostat}'.")

    with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tarball:
        tmp_name = tarball.name

    logger.info(f"Using temporary file: {tmp_name}")
    progress_retrieve(url_eurostat, tmp_name, disable=disable_progress)

    logger.info("Extracting Eurostat data.")
    with zipfile.ZipFile(tmp_name, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    logger.info(f"Eurostat data available in '{to_fn}'.")

    Path(tmp_name).unlink(missing_ok=True)
