# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build historical annual ammonia production per country in ktonNH3/a.
"""

import country_converter as coco
import pandas as pd

cc = coco.CountryConverter()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_ammonia_production")

    ammonia = pd.read_excel(
        snakemake.input.usgs,
        sheet_name="T12",
        skiprows=5,
        header=0,
        index_col=0,
        skipfooter=19,
    )

    ammonia.index = cc.convert(ammonia.index, to="iso2")

    years = [str(i) for i in range(2013, 2018)]
    countries = ammonia.index.intersection(snakemake.config["countries"])
    ammonia = ammonia.loc[countries, years].astype(float)

    # convert from ktonN to ktonNH3
    ammonia *= 17 / 14

    ammonia.index.name = "ktonNH3/a"

    ammonia.to_csv(snakemake.output.ammonia_production)
