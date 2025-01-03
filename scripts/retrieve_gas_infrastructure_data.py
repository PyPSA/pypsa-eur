# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve gas infrastructure data from
https://zenodo.org/records/4767098/files/IGGIELGN.zip.
"""

import logging
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

        snakemake = mock_snakemake("retrieve_gas_network_data")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    url = "https://zenodo.org/records/4767098/files/IGGIELGN.zip"

    # Save locations
    zip_fn = Path(f"{rootpath}/IGGIELGN.zip")
    to_fn = Path(rootpath) / Path(snakemake.output[0]).parent.parent

    logger.info(f"Downloading databundle from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, zip_fn, disable=disable_progress)

    validate_checksum(zip_fn, url)

    logger.info("Extracting databundle.")
    zipfile.ZipFile(zip_fn).extractall(to_fn)

    zip_fn.unlink()

    logger.info(f"Gas infrastructure data available in '{to_fn}'.")
