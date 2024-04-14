# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Retrieve and extract eurostat energy balances data.
"""


import logging
import zipfile
import gzip
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
    url_eurostat = "https://ec.europa.eu/eurostat/documents/38154/4956218/Balances-December2022.zip/f7cf0d19-5c0f-60ad-4e48-098a5ddd6e48?t=1671184070589"
    tarball_fn = Path(f"{rootpath}/data/eurostat/eurostat_2023.zip")
    to_fn = Path(
        f"{rootpath}/data/eurostat/eurostat-energy_balances-april_2023_edition/"
    )

    logger.info(f"Downloading Eurostat data from '{url_eurostat}'.")
    progress_retrieve(url_eurostat, tarball_fn, disable=disable_progress)

    logger.info("Extracting Eurostat data.")
    with zipfile.ZipFile(tarball_fn, "r") as zip_ref:
        zip_ref.extractall(to_fn)

    logger.info(f"Eurostat data available in '{to_fn}'.")

    url_eurostat_household = "https://ec.europa.eu/eurostat/api/dissemination/sdmx/3.0/data/dataflow/ESTAT/nrg_d_hhq/1.0/*.*.*.*.*?c[freq]=A&c[nrg_bal]=FC_OTH_HH_E,FC_OTH_HH_E_SH,FC_OTH_HH_E_WH,FC_OTH_HH_E_CK&c[siec]=TOTAL&c[unit]=TJ&c[geo]=EU27_2020,EA20,BE,BG,CZ,DK,DE,EE,IE,EL,ES,FR,HR,IT,CY,LV,LT,LU,HU,MT,NL,AT,PL,PT,RO,SI,SK,FI,SE,NO,UK,BA,MD,MK,AL,RS,UA,XK,GE&compress=true&format=csvdata&formatVersion=2.0&c[time]=2021,2020,2019,2018,2017,2016,2015,2014,2013,2012,2011,2010"
    tarball_fn = Path(f"{rootpath}/data/eurostat/eurostat_household.gz")
    to_fn = Path(
        f"{rootpath}/data/eurostat/eurostat-household_energy_balances-february_2024.csv"
    )

    logger.info(f"Downloading Eurostats' disaggregated household energy balances data from '{url_eurostat_household}'.")
    progress_retrieve(url_eurostat_household, tarball_fn, disable=disable_progress)

    logger.info("Extracting Eurostat's disaggregated household energy balance data.")
    with gzip.open(tarball_fn, 'rb') as f_in, open(to_fn, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    logger.info(f"Eurostat's disaggregated household energy balance data available in '{to_fn}'.")
