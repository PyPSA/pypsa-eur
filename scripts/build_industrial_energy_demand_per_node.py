# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build industrial energy demand per model region.

Description
-------
This rule aggregates the energy demand of the industrial sectors per model region.
For each bus, the following carriers are considered:
- electricity
- coal
- coke
- solid biomass
- methane
- hydrogen
- low-temperature heat
- naphtha
- ammonia
- process emission
- process emission from feedstock

which can later be used as values for the industry load.
"""

import logging

import pandas as pd
from _helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_energy_demand_per_node",
            clusters=48,
            planning_horizons=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # import ratios
    fn = snakemake.input.industry_sector_ratios
    sector_ratios = pd.read_csv(fn, header=[0, 1], index_col=0)

    # material demand per node and industry (Mton/a)
    fn = snakemake.input.industrial_production_per_node
    nodal_production = pd.read_csv(fn, index_col=0) / 1e3

    # energy demand today to get current electricity
    fn = snakemake.input.industrial_energy_demand_per_node_today
    nodal_today = pd.read_csv(fn, index_col=0)

    nodal_sector_ratios = pd.concat(
        {node: sector_ratios[node[:2]] for node in nodal_production.index}, axis=1
    )

    nodal_production_stacked = nodal_production.stack()
    nodal_production_stacked.index.names = [None, None]

    # final energy consumption per node and industry (TWh/a)
    nodal_df = (
        (nodal_sector_ratios.multiply(nodal_production_stacked))
        .T.groupby(level=0)
        .sum()
    )

    rename_sectors = {
        "elec": "electricity",
        "biomass": "solid biomass",
        "heat": "low-temperature heat",
    }
    nodal_df.rename(columns=rename_sectors, inplace=True)

    nodal_df["current electricity"] = nodal_today["electricity"]

    nodal_df.index.name = "TWh/a (MtCO2/a)"

    fn = snakemake.output.industrial_energy_demand_per_node
    nodal_df.to_csv(fn, float_format="%.2f")
