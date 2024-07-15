# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve monthly fuel prices from Destatis.
"""

import logging
from pathlib import Path

from _helpers import configure_logging, retrieve_file, set_scenario_config

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_gdp_uamd")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

dict_urls = dict(
    {
        "gdp": "https://datadryad.org/stash/downloads/file_stream/241947",
        "ppp": "https://github.com/ecohealthalliance/sars_cov_risk/releases/download/v2.0.1/ppp_2020_1km_Aggregated.tif",
    }
)

# Download and validate each dataset
for key, path in snakemake.output.items():
    retrieve_file(dict_urls[key], path)
