# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build industrial energy demand per model region.
"""

import pandas as pd

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_energy_demand_per_node",
            weather_year="",
            simpl="",
            clusters=48,
            planning_horizons=2030,
        )

    # import EU ratios df as csv
    fn = snakemake.input.industry_sector_ratios
    industry_sector_ratios = pd.read_csv(fn, index_col=0)

    # material demand per node and industry (kton/a)
    fn = snakemake.input.industrial_production_per_node
    nodal_production = pd.read_csv(fn, index_col=0)

    # energy demand today to get current electricity
    fn = snakemake.input.industrial_energy_demand_per_node_today
    nodal_today = pd.read_csv(fn, index_col=0)

    # final energy consumption per node and industry (TWh/a)
    nodal_df = nodal_production.dot(industry_sector_ratios.T)

    # convert GWh to TWh and ktCO2 to MtCO2
    nodal_df *= 0.001

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
