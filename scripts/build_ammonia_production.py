# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build historical annual ammonia production per country in ktonNH3/a.

Inputs
-------

- ``data/bundle-sector/myb1-2017-nitro.xls``

Outputs
-------

- ``resources/ammonia_production.csv``

Description
-------

This functions takes data from the `Minerals Yearbook <https://www.usgs.gov/centers/national-minerals-information-center/nitrogen-statistics-and-information>`_ (June 2024) published by the US Geological Survey (USGS) and the National Minerals Information Center and extracts the annual ammonia production per country in ktonN/a. The data is converted to ktonNH3/a.
"""

import country_converter as coco
import pandas as pd
from _helpers import set_scenario_config

cc = coco.CountryConverter()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ammonia_production")

    set_scenario_config(snakemake)

    ammonia = pd.read_excel(
        snakemake.input.usgs,
        sheet_name="T12",
        skiprows=5,
        header=0,
        index_col=0,
        skipfooter=19,
        na_values=["--"],
    )

    ammonia.index = cc.convert(ammonia.index, to="iso2")

    years = [str(i) for i in range(2013, 2018)]

    ammonia = ammonia[years]

    # convert from ktonN to ktonNH3
    ammonia *= 17 / 14

    ammonia.index.name = "ktonNH3/a"

    ammonia.to_csv(snakemake.output.ammonia_production)
