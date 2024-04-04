# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build specific energy consumption by carrier and industries and by country,
that interpolates between the current average energy consumption (from
2015-2020) and the ideal future best-in-class consumption.
"""

import pandas as pd
from prepare_sector_network import get


def build_industry_sector_ratios_intermediate():

    # in TWh/a
    demand = pd.read_csv(
        snakemake.input.industrial_energy_demand_per_country_today,
        header=[0, 1],
        index_col=0,
    )

    # in Mt/a
    production = (
        pd.read_csv(snakemake.input.industrial_production_per_country, index_col=0)
        / 1e3
    ).stack()
    production.index.names = [None, None]

    # in MWh/t
    future_sector_ratios = pd.read_csv(
        snakemake.input.industry_sector_ratios, index_col=0
    )

    today_sector_ratios = demand.div(production, axis=1)

    today_sector_ratios.dropna(how="all", axis=1, inplace=True)

    rename = {
        "waste": "biomass",
        "electricity": "elec",
        "solid": "coke",
        "gas": "methane",
        "other": "biomass",
        "liquid": "naphtha",
    }
    today_sector_ratios = today_sector_ratios.rename(rename).groupby(level=0).sum()

    fraction_future = get(params["sector_ratios_fraction_future"], year)

    intermediate_sector_ratios = {}
    for ct, group in today_sector_ratios.T.groupby(level=0):
        today_sector_ratios_ct = (
            group.droplevel(0)
            .T.reindex_like(future_sector_ratios)
        )
        missing_mask = (today_sector_ratios_ct.isna().all())
        today_sector_ratios_ct.loc[:, missing_mask] = future_sector_ratios.loc[:, missing_mask]
        today_sector_ratios_ct.loc[:, ~missing_mask] = today_sector_ratios_ct.loc[:, ~missing_mask].fillna(0)
        intermediate_sector_ratios[ct] = (
            today_sector_ratios_ct * (1 - fraction_future)
            + future_sector_ratios * fraction_future
        )
    intermediate_sector_ratios = pd.concat(intermediate_sector_ratios, axis=1)

    intermediate_sector_ratios.to_csv(snakemake.output.industry_sector_ratios)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_sector_ratios_intermediate",
            planning_horizons="2030",
        )

    year = int(snakemake.wildcards.planning_horizons[-4:])

    params = snakemake.params.industry

    build_industry_sector_ratios_intermediate()
