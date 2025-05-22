# SPDX-FileCopyrightText: 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.14230568.svg
  :target: https://doi.org/10.5281/zenodo.14230568

The data bundle contains input data for the 2024 TYNDP scenario building process.

This rule downloads the TYNDP data bundle from `zenodo
<https://doi.org/10.5281/zenodo.14230568>` and extracts it in the ``data``
subdirectory, such that all files of the bundle are stored in the
``data/tyndp_2024_bundle`` subdirectory.

**Outputs**

- ``data/tyndp_2024_bundle``: input data for TYNDP 2024 scenario building

"""

import logging
import os
import zipfile

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

# Define the base URL
url = "https://zenodo.org/records/14230568/files/TYNDP_2024_data_bundle.zip"

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_tyndp_bundle")
        rootpath = ".."
    else:
        rootpath = "."

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    to_fn = snakemake.output.dir
    to_fn_zp = to_fn + ".zip"

    # download .zip file
    logger.info(f"Downloading TYNDP data bundle from '{url}'.")
    progress_retrieve(url, to_fn_zp, disable=disable_progress)

    # extract
    logger.info("Extracting TYNDP data bundle.")
    with zipfile.ZipFile(to_fn_zp, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    # remove .zip file
    os.remove(to_fn_zp)

    logger.info(f"TYNDP data bundle available in '{to_fn}'.")
