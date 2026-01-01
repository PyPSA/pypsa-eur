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

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)

sectors_copied_from_exogenous = [
    # also mainly has a heat demand, but can be electrified by methods unavailable in other processes.
    # Therefore, here exogenously electrified as long as individual sectors are not disaggregated.
    "HVC (chemical recycling)",
    "HVC (mechanical recycling)",
    "DRI + Electric arc",
]

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_energy_demand_per_node",
            clusters=48,
            planning_horizons=2030,
        )
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # import exogenous ratios
    fn = snakemake.input.industry_sector_ratios
    sector_ratios_exogenous = pd.read_csv(fn, header=[0, 1], index_col=0)

    # import endogenous ratios
    fn = snakemake.input.industry_sector_ratios_endogenous
    sector_ratios_endogenous = pd.read_csv(fn, header=[0, 1], index_col=0)

    cols = sector_ratios_exogenous.columns[
        sector_ratios_exogenous.columns.get_level_values(1).isin(
            sectors_copied_from_exogenous
        )
    ]
    sector_ratios_endogenous = pd.concat(
        [sector_ratios_endogenous, sector_ratios_exogenous.loc[:, cols]], axis=1
    )

    # material demand per node and industry (Mton/a)
    fn = snakemake.input.industrial_production_per_node
    nodal_production = pd.read_csv(fn, index_col=0) / 1e3

    # energy demand today to get current electricity
    fn = snakemake.input.industrial_energy_demand_per_node_today
    nodal_today = pd.read_csv(fn, index_col=0)

    nodal_sector_ratios_exogenous = pd.concat(
        {node: sector_ratios_exogenous[node[:2]] for node in nodal_production.index},
        axis=1,
    )

    nodal_sector_ratios_endogenous = pd.concat(
        {node: sector_ratios_endogenous[node[:2]] for node in nodal_production.index},
        axis=1,
    )

    nodal_production_stacked = nodal_production.stack()
    nodal_production_stacked.index.names = [None, None]

    # final energy consumption per node and industry (TWh/a)
    nodal_df_exogenous = (
        (nodal_sector_ratios_exogenous.multiply(nodal_production_stacked))
        .T.groupby(level=0)
        .sum()
    )

    nodal_df_endogenous = (
        (nodal_sector_ratios_endogenous.multiply(nodal_production_stacked))
        .T.groupby(level=0)
        .sum()
    )

    rename_sectors = pd.Series(
        {
            "elec": "electricity",
            "biomass": "solid biomass",
        }
    )

    nodal_df_exogenous.rename(columns=rename_sectors, inplace=True)
    nodal_df_endogenous.rename(columns=rename_sectors, inplace=True)

    nodal_df_exogenous.rename(
        columns={
            "heat": "low-temperature heat",
        },
        inplace=True,
    )

    grouper = {
        "heat<100": "low-temperature heat",
        "heat": "low-temperature heat",
    }

    def partial_group(df, grouper):
        return pd.concat([df.groupby(grouper).sum(), df.drop(grouper)])

    nodal_df_endogenous = partial_group(nodal_df_endogenous.T, grouper).T

    logger.warning(
        "what about methanol, process emissions from feedstock and current electricity?"
    )
    nodal_df_exogenous["current electricity"] = nodal_today["electricity"]
    nodal_df_endogenous["current electricity"] = nodal_today["electricity"]

    nodal_df = pd.concat(
        [nodal_df_exogenous, nodal_df_endogenous],
        axis=1,
        keys=["exogenous", "endogenous"],
    )

    idx = pd.IndexSlice
    nodal_df.index.name = "TWh/a (MtCO2/a)"

    fn = snakemake.output.industrial_energy_demand_per_node
    nodal_df.to_csv(fn, float_format="%.2f")
