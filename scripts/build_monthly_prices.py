# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script extracts monthly fuel prices of oil, gas, coal and lignite, as well
as CO2 prices.

Description
-----------

The rule :mod:`build_monthly_prices` collects monthly fuel prices and CO2 prices
and translates them from different input sources to pypsa syntax

Data sources:
    [1] Fuel price index. Destatis
    https://www.destatis.de/EN/Home/_node.html
    [2] average annual fuel price lignite, ENTSO-E
    https://2020.entsos-tyndp-scenarios.eu/fuel-commodities-and-carbon-prices/
    [3] CO2 Prices, Emission spot primary auction, EEX
    https://www.eex.com/en/market-data/environmental-markets/eua-primary-auction-spot-download


Data was accessed at 16.5.2023
"""

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


# keywords in datasheet
keywords = {
    "coal": " GP09-051 Hard coal",
    "lignite": " GP09-052 Lignite and lignite briquettes",
    "oil": " GP09-0610 10 Mineral oil, crude",
    "gas": "GP09-062 Natural gas",
}

# sheet names to pypsa syntax
sheet_name_map = {
    "coal": "5.1 Hard coal and lignite",
    "lignite": "5.1 Hard coal and lignite",
    "oil": "5.2 Mineral oil",
    "gas": "5.3.1 Natural gas - indices",
}


# import fuel price 2015 in Eur/MWh
# source lignite, price for 2020, scaled by price index, ENTSO-E [3]
price_2020 = (
    pd.Series({"coal": 3.0, "oil": 10.6, "gas": 5.6, "lignite": 1.1}) * 3.6
)  # Eur/MWh

# manual adjustment of coal price
price_2020["coal"] = 2.4 * 3.6
price_2020["lignite"] = 1.6 * 3.6


def get_fuel_price():
    price = {}
    for carrier, keyword in keywords.items():
        sheet_name = sheet_name_map[carrier]
        df = pd.read_excel(
            snakemake.input.fuel_price_raw,
            sheet_name=sheet_name,
            index_col=0,
            skiprows=6,
            nrows=18,
        )
        df = df.dropna(axis=0).iloc[:, :12]
        start, end = df.index[0], str(int(df.index[-1][:4]) + 1)
        df = df.stack()
        df.index = pd.date_range(start=start, end=end, freq="MS", inclusive="left")
        scale = price_2020[carrier] / df["2020"].mean()  # scale to 2020 price
        df = df.mul(scale)
        price[carrier] = df

    return pd.concat(price, axis=1)


def get_co2_price():
    # emission price
    co2_price = pd.read_excel(snakemake.input.co2_price_raw, index_col=1, header=5)
    return co2_price["Auction Price €/tCO2"]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_monthly_prices")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    fuel_price = get_fuel_price()
    fuel_price.to_csv(snakemake.output.fuel_price)

    co2_price = get_co2_price()
    co2_price.to_csv(snakemake.output.co2_price)
