# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve cost data from ``technology-data``.
"""

import logging
from pathlib import Path

from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_cost_data", year=2030)
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    version = snakemake.params.version
    if "/" in version:
        baseurl = f"https://raw.githubusercontent.com/{version}/outputs"
    else:
        baseurl = f"https://raw.githubusercontent.com/PyPSA/technology-data/{version}/outputs/"
    filepath = Path(snakemake.output[0])
    url = baseurl + filepath.name

    print(url)

    to_fn = Path(rootpath) / filepath

    print(to_fn)

    logger.info(f"Downloading technology data from '{url}'.")
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, to_fn, disable=disable_progress)

    logger.info(f"Technology data available at at {to_fn}")
