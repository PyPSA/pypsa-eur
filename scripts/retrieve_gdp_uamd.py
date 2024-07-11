# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve monthly fuel prices from Destatis.
"""

import logging
from pathlib import Path

from _helpers import retrieve_file

logger = logging.getLogger(__name__)


def retrieve(url, destination):

    logger.info(f"Downloading file from '{url}'.")
    retrieve_file(url, destination)
    logger.info("File downloaded and validated.")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_gdp_uamd")
        rootpath = ".."
    else:
        rootpath = "."

datasets = [
    # GDP_PPP_30arcsec_v3.nc: raw dataset. Available at: [M. Kummu, M. Taka, J. H. A. Guillaume. (2020), Data from: Gridded global datasets for Gross Domestic Product and Human Development Index over 1990-2015, Dryad, Dataset. doi: https://doi.org/10.5061/dryad.dk1j0]
    (
        "https://datadryad.org/stash/downloads/file_stream/241947",
        "GDP_per_capita_PPP_1990_2015_v2.nc",
    ),
    # ppp_2020_1km_Aggregated.tif: raw dataset. Available at: https://data.humdata.org/dataset/
    (
        "https://data.worldpop.org/GIS/Population/Global_2000_2020/2020/0_Mosaicked/ppp_2020_1km_Aggregated.tif",
        "ppp_2020_1km_Aggregated.tif",
    ),
]

# Download and validate each dataset
for url, filename in datasets:
    file_path = rootpath / "data" / filename
    retrieve(url, file_path)
