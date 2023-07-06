# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 10:37:35 2023.

This script extracts monthly fuel prices of oil, gas, coal and lignite,
as well as CO2 prices


Inputs
------
- ``data/energy-price-trends-xlsx-5619002.xlsx``: energy price index of fossil fuels
- ``emission-spot-primary-market-auction-report-2019-data.xls``: CO2 Prices spot primary auction


Outputs
-------

- ``data/validation/monthly_fuel_price.csv``
- ``data/validation/CO2_price_2019.csv``


Description
-----------

The rule :mod:`build_monthly_prices` collects monthly fuel prices and CO2 prices
and translates them from different input sources to pypsa syntax

Data sources:
    [1] Fuel price index. Destatis
    https://www.destatis.de/EN/Home/_node.html
    [2] average annual import price (coal, gas, oil) Agora, slide 24
    https://static.agora-energiewende.de/fileadmin/Projekte/2019/Jahresauswertung_2019/A-EW_German-Power-Market-2019_Summary_EN.pdf
    [3] average annual fuel price lignite, ENTSO-E
    https://2020.entsos-tyndp-scenarios.eu/fuel-commodities-and-carbon-prices/
    [4] CO2 Prices, Emission spot primary auction, EEX
    https://www.eex.com/en/market-data/environmental-markets/eua-primary-auction-spot-download


Data was accessed at 16.5.2023
"""

import logging

import pandas as pd
from _helpers import configure_logging

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
# source for coal, oil, gas, Agora, slide 24 [2]
# source lignite, price for 2020, scaled by price index, ENTSO-E [3]
price_2015 = {"coal": 8.3, "oil": 30.6, "gas": 20.6, "lignite": 3.8}  # 2020 3.96/1.04


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
        df = df.mul(price_2015[carrier] / 100)
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

    fuel_price = get_fuel_price()
    fuel_price.to_csv(snakemake.output.fuel_price)

    co2_price = get_co2_price()
    co2_price.to_csv(snakemake.output.co2_price)
