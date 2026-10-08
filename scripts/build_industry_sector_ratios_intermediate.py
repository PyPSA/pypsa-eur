# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build specific energy consumption per industrial subsector, carrier and country for a planning horizon in MWh/t.

Today's ratios are derived per country by dividing the energy demand from [build_industrial_energy_demand_per_country_today][] by the production from [build_industrial_production_per_country][], with the carriers mapped onto those of the best-case ratios (solid to coke, gas to methane, liquid to naphtha, waste and other to biomass). These are then interpolated linearly towards the best-case ratios from [build_industry_sector_ratios][] using `industry: sector_ratios_fraction_future` for the horizon year: a fraction of 0 keeps today's ratios, 1 uses the best-case ratios. Subsectors that do not exist today, such as DRI + Electric arc or HVC recycling, and any missing country or carrier values take the best-case ratios.
"""

import logging

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config
from scripts.prepare_sector_network import get

logger = logging.getLogger(__name__)


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

    today_sector_ratios = demand.div(production, axis=1).replace([np.inf, -np.inf], 0)

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
        today_sector_ratios_ct = group.droplevel(0).T.reindex_like(future_sector_ratios)
        missing_mask = today_sector_ratios_ct.isna().all()
        today_sector_ratios_ct.loc[:, missing_mask] = future_sector_ratios.loc[
            :, missing_mask
        ]
        today_sector_ratios_ct.loc[:, ~missing_mask] = today_sector_ratios_ct.loc[
            :, ~missing_mask
        ].fillna(future_sector_ratios)
        intermediate_sector_ratios[ct] = (
            today_sector_ratios_ct * (1 - fraction_future)
            + future_sector_ratios * fraction_future
        )

    intermediate_sector_ratios = pd.concat(intermediate_sector_ratios, axis=1)

    intermediate_sector_ratios.to_csv(snakemake.output.industry_sector_ratios)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_sector_ratios_intermediate",
            horizon="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    year = int(snakemake.wildcards.horizon)

    params = snakemake.params.industry

    build_industry_sector_ratios_intermediate()
