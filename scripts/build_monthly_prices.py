# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script extracts monthly fuel prices of oil, gas and coal.

Description
-----------

The rule :mod:`build_monthly_prices` collects monthly fuel prices
and translates them from different input sources to pypsa syntax.
"""

import logging

import pandas as pd
from pydeflate import imf_gdp_deflate, set_pydeflate_path

from scripts._helpers import configure_logging, set_scenario_config

set_pydeflate_path("../data/pydeflate/")

logger = logging.getLogger(__name__)

MMBTU_PER_MWH = 3.41214
BBL_PER_MWH = 0.5883
METRIC_TON_PER_MWH_COAL = 0.1433

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_fossil_fuel_prices")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    df = pd.read_excel(
        snakemake.input.fuel_price_raw,
        skiprows=[0, 1, 2, 3, 5],
        header=0,
        index_col=0,
        sheet_name="Monthly Prices",
        parse_dates=True,
        date_format="%YM%m",
        na_values=["…"],
    )

    COLUMNS = {
        "Crude oil, Brent": "oil",
        "Coal, South African **": "coal",
        "Natural gas, Europe": "gas",
    }
    df = df.rename(columns=COLUMNS)[list(COLUMNS.values())]

    # Convert from nominal to real prices (deflate to base year 2020)
    df["year"] = df.index.year
    df["iso_code"] = "DEU"
    for col in COLUMNS.values():
        df[col] = imf_gdp_deflate(
            df,
            value_column=col,
            source_currency="USD",
            target_currency="EUR",
            base_year=2020,
        )["value"].values
    df = df[df.index.year >= 1999]  # only available from 1999 onwards
    df = df.drop(columns=["year", "iso_code"])

    df["oil"] *= BBL_PER_MWH
    df["gas"] *= MMBTU_PER_MWH
    df["coal"] *= METRIC_TON_PER_MWH_COAL

    # rolling mean for smoothing
    window = snakemake.params.rolling_window
    df = df.rolling(window=window, center=True, min_periods=1).mean()

    df.to_csv(snakemake.output.fuel_price)
