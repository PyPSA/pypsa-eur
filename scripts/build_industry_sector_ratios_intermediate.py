# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build specific energy consumption by carrier and industries and by country,
that interpolates between the current average energy consumption (from 2015-2020)
and the ideal future best-in-class consumption.
"""

import pandas as pd

from prepare_sector_network import get

def build_industry_sector_ratios_intermediate():

    # in TWh/a
    demand = pd.read_csv(snakemake.input.industrial_energy_demand_per_country_today,
                         header=[0,1],
                         index_col=0)

    # in Mt/a
    production = pd.read_csv(snakemake.input.industrial_production_per_country,
                             index_col=0) / 1e3
    production = production.unstack().swaplevel()

    # in MWh/t
    future_sector_ratios = pd.read_csv(snakemake.input.industry_sector_ratios,
                                       index_col=0)

    production.index.names = [None,None]

    today_sector_ratios = demand.div(production, axis=1)

    today_sector_ratios.drop(columns=today_sector_ratios.columns[today_sector_ratios.isna().all()],
                             inplace=True)

    rename = pd.Series(today_sector_ratios.index,
                       today_sector_ratios.index)
    rename["waste"] = "biomass"
    rename["electricity"] = "elec"
    rename["solid"] = "coke"
    rename["gas"] = "methane"
    rename["other"] = "biomass"
    rename["liquid"] = "naphtha"

    today_sector_ratios.rename(rename,
                               inplace=True)


    fraction_future = get(params["sector_ratios_fraction_future"], year)

    intermediate_sector_ratios = {}

    for ct in today_sector_ratios.columns.unique(level=0):

        intermediate_sector_ratio = future_sector_ratios.copy()

        intermediate_sector_ratio.loc[today_sector_ratios[ct].index,today_sector_ratios[ct].columns] = (fraction_future*intermediate_sector_ratio.loc[today_sector_ratios[ct].index,today_sector_ratios[ct].columns]
                                                                                                        + (1 - fraction_future)*today_sector_ratios[ct])
        intermediate_sector_ratios[ct] = intermediate_sector_ratio

    intermediate_sector_ratios = pd.concat(intermediate_sector_ratios, axis=1)

    intermediate_sector_ratios.to_csv(snakemake.output.industry_sector_ratios)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industry_sector_ratios_intermediate")

    year = int(snakemake.wildcards.planning_horizons[-4:])

    params = snakemake.params.industry

    build_industry_sector_ratios_intermediate()
