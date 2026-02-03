# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script extracts monthly fuel prices of oil, gas, coal and lignite.

Description
-----------

The rule :mod:`build_monthly_prices` collects monthly fuel prices
and translates them from different input sources to pypsa syntax.
"""

import logging

import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

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
    )
