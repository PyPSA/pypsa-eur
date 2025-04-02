# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3517934.svg
   :target: https://doi.org/10.5281/zenodo.3517934

The data bundle contains common GIS datasets like NUTS3 shapes, EEZ shapes,
CORINE Landcover, Natura 2000 and also electricity specific summary statistics
like historic per country yearly totals of hydro generation, GDP and population
data on NUTS3 levels and energy balances.

This rule downloads the data bundle from `zenodo
<https://doi.org/10.5281/zenodo.3517934>`_ and extracts it in the ``data``
sub-directory, such that all files of the bundle are stored in the
``data/bundle`` subdirectory.
"""

import logging
import tarfile
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

    url = "https://zenodo.org/records/14732152/files/bundle.tar.xz"

    tarball_fn = Path(f"{rootpath}/bundle.tar.xz")
    to_fn = Path(rootpath) / Path(snakemake.output[0]).parent.parent

    logger.info(f"Downloading databundle from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, tarball_fn, disable=disable_progress)

    validate_checksum(tarball_fn, url)

    logger.info("Extracting databundle.")
    tarfile.open(tarball_fn).extractall(to_fn)

    logger.info("Unlinking tarball.")
    tarball_fn.unlink()

    logger.info(f"Databundle available in '{to_fn}'.")
