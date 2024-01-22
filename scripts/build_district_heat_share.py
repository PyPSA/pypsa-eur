# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build district heat shares at each node, depending on investment year.
"""

import logging

import pandas as pd
from prepare_sector_network import get

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_district_heat_share",
            simpl="",
            clusters=48,
            planning_horizons="2050",
        )

    investment_year = int(snakemake.wildcards.planning_horizons[-4:])

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    district_heat_share = pd.read_csv(snakemake.input.district_heat_share, index_col=0)[
        "district heat share"
    ]

    # make ct-based share nodal
    district_heat_share = district_heat_share.loc[pop_layout.ct]
    district_heat_share.index = pop_layout.index

    # total urban population per country
    ct_urban = pop_layout.urban.groupby(pop_layout.ct).sum()

    # distribution of urban population within a country
    pop_layout["urban_ct_fraction"] = pop_layout.urban / pop_layout.ct.map(ct_urban.get)

    # fraction of node that is urban
    urban_fraction = pop_layout.urban / pop_layout[["rural", "urban"]].sum(axis=1)

    # maximum potential of urban demand covered by district heating
    central_fraction = snakemake.config["sector"]["district_heating"]["potential"]

    # district heating share at each node
    dist_fraction_node = (
        district_heat_share * pop_layout["urban_ct_fraction"] / pop_layout["fraction"]
    )

    # if district heating share larger than urban fraction -> set urban
    # fraction to district heating share
    urban_fraction = pd.concat([urban_fraction, dist_fraction_node], axis=1).max(axis=1)

    # difference of max potential and today's share of district heating
    diff = (urban_fraction * central_fraction) - dist_fraction_node
    progress = get(
        snakemake.config["sector"]["district_heating"]["progress"], investment_year
    )
    dist_fraction_node += diff * progress
    logger.info(
        f"Increase district heating share by a progress factor of {progress:.2%} "
        f"resulting in new average share of {dist_fraction_node.mean():.2%}"
    )

    df = pd.DataFrame(
        {
            "original district heat share": district_heat_share,
            "district fraction of node": dist_fraction_node,
            "urban fraction": urban_fraction,
        },
        dtype=float,
    )

    df.to_csv(snakemake.output.district_heat_share)
