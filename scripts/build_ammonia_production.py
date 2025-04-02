# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build historical annual ammonia production per country in ktonNH3/a.

Description
-------

This functions takes data from the `Minerals Yearbook <https://www.usgs.gov/centers/national-minerals-information-center/nitrogen-statistics-and-information>`_
 (July 2024) published by the US Geological Survey (USGS) and the National Minerals Information Center and extracts the annual ammonia production per country in ktonN/a. The data is converted to ktonNH3/a.
"""

import logging

import country_converter as coco
import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ammonia_production")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    ammonia = pd.read_excel(
        snakemake.input.usgs,
        sheet_name="T12",
        skiprows=5,
        header=0,
        index_col=0,
        skipfooter=7,
        na_values=["--"],
    )

    ammonia.index = cc.convert(ammonia.index, to="iso2")

    years = [str(i) for i in range(2018, 2023)]

    ammonia = ammonia.rename(columns=lambda x: str(x))[years]

    # convert from ktonN to ktonNH3
    ammonia *= 17 / 14

    ammonia.index.name = "ktonNH3/a"

    ammonia.to_csv(snakemake.output.ammonia_production)
