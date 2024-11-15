# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
"""

import logging
from pathlib import Path

import requests
from _helpers import configure_logging, progress_retrieve, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_heat_source_utilisation_potentials")
        rootpath = ".."
    else:
        rootpath = "."
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # license: https://creativecommons.org/licenses/by/4.0/
    # download the data in url
    heat_source = snakemake.params["heat_source"]
    filepath = Path(snakemake.output[0])
    if not filepath.parent.exists():
        filepath.parent.mkdir(parents=True)

    url = f"https://fordatis.fraunhofer.de/bitstream/fordatis/341.3/10/{snakemake.params.fraunhofer_heat_utilisation_potentials[heat_source]['key']}.gpkg"

    logger.info(
        f"Downloading heat source utilisation potential data for {heat_source} from '{url}'."
    )
    disable_progress = snakemake.config["run"].get("disable_progressbar", False)
    progress_retrieve(url, filepath, disable=disable_progress)

    logger.info(f"Data available at at {filepath}")
