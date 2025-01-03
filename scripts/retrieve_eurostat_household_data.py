# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract eurostat household energy balances data.
"""

import gzip
import logging
import shutil
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

    url_eurostat_household = "https://ec.europa.eu/eurostat/databrowser-backend/api/extraction/1.0/LIVE/false/sdmx/csv/nrg_d_hhq__custom_11480365?startPeriod=2013&endPeriod=2022&i&compressed=true"
    tarball_fn = Path(f"{rootpath}/data/eurostat/eurostat_household.gz")
    to_fn = Path(
        f"{rootpath}/data/eurostat/eurostat-household_energy_balances-february_2024.csv"
    )

    logger.info(
        f"Downloading Eurostats' disaggregated household energy balances data from '{url_eurostat_household}'."
    )
    progress_retrieve(url_eurostat_household, tarball_fn, disable=disable_progress)

    logger.info("Extracting Eurostat's disaggregated household energy balance data.")
    with gzip.open(tarball_fn, "rb") as f_in, open(to_fn, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    logger.info(
        f"Eurostat's disaggregated household energy balance data available in '{to_fn}'."
    )
