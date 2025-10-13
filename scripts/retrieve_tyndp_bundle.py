# SPDX-FileCopyrightText: 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
The retrieved data includes the node list and the reference grids from the Ten-Year Network Development Plan (TYNDP) 2024.

This rule downloads the TYNDP data bundle from the `ENTSOs website
<https://2024.entsos-tyndp-scenarios.eu/download/>` and extracts it in the ``data/tyndp_2024_bundle``
subdirectory.

**Outputs**

- ``data/tyndp_2024_bundle/Line data/ReferenceGrid_Electricity.xlsx``: reference grid from TYNDP 2024
- ``data/tyndp_2024_bundle/Nodes/LIST OF NODES.xlsx``: list of nodes from TYNDP 2024

"""

import logging
import os
import shutil
import zipfile
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Define the base URL
url_electricity = (
    "https://2024-data.entsos-tyndp-scenarios.eu/files/scenarios-inputs/Line-data.zip"
)
url_buses = (
    "https://2024-data.entsos-tyndp-scenarios.eu/files/scenarios-inputs/Nodes.zip"
)


def retrieve_bundle(url: str, to_fn: str, disable_progress: bool = False):
    to_fn_zp = Path(to_fn, Path(url).name)

    # download .zip file
    logger.info(f"Downloading TYNDP data bundle from '{url}'.")
    progress_retrieve(url, to_fn_zp, disable=disable_progress)

    # extract
    logger.info("Extracting TYNDP data bundle.")
    with zipfile.ZipFile(to_fn_zp, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    # remove .zip file and __MACOSX
    os.remove(to_fn_zp)
    shutil.rmtree(Path(to_fn, "__MACOSX"), ignore_errors=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_tyndp_bundle")

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)

    to_fn = Path(snakemake.output.reference_grid).parents[1]
    urls = [
        url_electricity,
        url_buses,
    ]

    # Retrieve TYNDP data
    tqdm_kwargs = {
        "ascii": False,
        "unit": " bundle",
        "total": len(urls),
        "desc": "Retrieving TYNDP data bundle",
        "disable": disable_progress,
    }

    for url in tqdm(urls, **tqdm_kwargs):
        retrieve_bundle(url, to_fn, disable_progress)

    logger.info(f"TYNDP data bundle available in '{to_fn}'.")
