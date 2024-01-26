# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Builds table of existing heat generation capacities for initial planning
horizon.
"""
import country_converter as coco
import numpy as np
import pandas as pd

cc = coco.CountryConverter()


def build_existing_heating():
    # retrieve existing heating capacities

    existing_heating = pd.read_csv(
        snakemake.input.existing_heating, index_col=0, header=0
    )

    # data for Albania, Montenegro and Macedonia not included in database
    existing_heating.loc["Albania"] = np.nan
    existing_heating.loc["Montenegro"] = np.nan
    existing_heating.loc["Macedonia"] = np.nan

    existing_heating.fillna(0.0, inplace=True)

    # convert GW to MW
    existing_heating *= 1e3

    existing_heating.index = cc.convert(existing_heating.index, to="iso2")

    # coal and oil boilers are assimilated to oil boilers
    existing_heating["oil boiler"] = (
        existing_heating["oil boiler"] + existing_heating["coal boiler"]
    )
    existing_heating.drop(["coal boiler"], axis=1, inplace=True)

    # distribute technologies to nodes by population
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    nodal_heating = existing_heating.loc[pop_layout.ct]
    nodal_heating.index = pop_layout.index
    nodal_heating = nodal_heating.multiply(pop_layout.fraction, axis=0)

    district_heat_info = pd.read_csv(snakemake.input.district_heat_share, index_col=0)
    dist_fraction = district_heat_info["district fraction of node"]
    urban_fraction = district_heat_info["urban fraction"]

    energy_layout = pd.read_csv(
        snakemake.input.clustered_pop_energy_layout, index_col=0
    )

    uses = ["space", "water"]
    sectors = ["residential", "services"]

    nodal_sectoral_totals = pd.DataFrame(dtype=float)

    for sector in sectors:
        nodal_sectoral_totals[sector] = energy_layout[
            [f"total {sector} {use}" for use in uses]
        ].sum(axis=1)

    nodal_sectoral_fraction = nodal_sectoral_totals.div(
        nodal_sectoral_totals.sum(axis=1), axis=0
    )

    nodal_heat_name_fraction = pd.DataFrame(dtype=float)

    nodal_heat_name_fraction["urban central"] = dist_fraction

    for sector in sectors:
        nodal_heat_name_fraction[f"{sector} rural"] = nodal_sectoral_fraction[
            sector
        ] * (1 - urban_fraction)
        nodal_heat_name_fraction[f"{sector} urban decentral"] = nodal_sectoral_fraction[
            sector
        ] * (urban_fraction - dist_fraction)

    nodal_heat_name_tech = pd.concat(
        {
            name: nodal_heating.multiply(nodal_heat_name_fraction[name], axis=0)
            for name in nodal_heat_name_fraction.columns
        },
        axis=1,
        names=["heat name", "technology"],
    )

    # move all ground HPs to rural, all air to urban

    for sector in sectors:
        nodal_heat_name_tech[(f"{sector} rural", "ground heat pump")] += (
            nodal_heat_name_tech[("urban central", "ground heat pump")]
            * nodal_sectoral_fraction[sector]
            + nodal_heat_name_tech[(f"{sector} urban decentral", "ground heat pump")]
        )
        nodal_heat_name_tech[(f"{sector} urban decentral", "ground heat pump")] = 0.0

        nodal_heat_name_tech[
            (f"{sector} urban decentral", "air heat pump")
        ] += nodal_heat_name_tech[(f"{sector} rural", "air heat pump")]
        nodal_heat_name_tech[(f"{sector} rural", "air heat pump")] = 0.0

    nodal_heat_name_tech[("urban central", "ground heat pump")] = 0.0

    nodal_heat_name_tech.to_csv(snakemake.output.existing_heating_distribution)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_existing_heating_distribution",
            simpl="",
            clusters=48,
            planning_horizons=2050,
        )

    build_existing_heating()
