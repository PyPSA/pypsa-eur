# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build specific energy consumption by carrier and industries and by country,
that interpolates between the current average energy consumption (from
2015-2020) and the ideal future best-in-class consumption.

Description
-------

The config["industry"]["sector_ratios_fraction_future"] parameter determines the progress towards the future best-in-class consumption.
For each bus, the following industry subcategories

- Electric arc
- DRI + Electric arc
- Integrated steelworks
- HVC
- HVC (mechanical recycling)
- HVC (chemical recycling)
- Ammonia
- Chlorine
- Methanol
- Other chemicals
- Pharmaceutical products etc.
- Cement
- Ceramics & other NMM
- Glass production
- Pulp production
- Paper production
- Printing and media reproduction
- Food, beverages and tobacco
- Alumina production
- Aluminium - primary production
- Aluminium - secondary production
- Other non-ferrous metals
- Transport equipment
- Machinery equipment
- Textiles and leather
- Wood and wood products
- Other Industrial Sectors

with the following carriers are considered:

- elec
- coal
- coke
- biomass
- methane
- hydrogen
- heat
- naphtha
- process emission
- process emission from feedstock
- (ammonia)

Unit of the output file is MWh/t.
"""

import logging

import numpy as np
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from prepare_sector_network import get

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
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_sector_ratios_intermediate",
            planning_horizons="2030",
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    year = int(snakemake.wildcards.planning_horizons)

    params = snakemake.params.industry

    build_industry_sector_ratios_intermediate()
