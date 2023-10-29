# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot average market prices and market values on map.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from _helpers import ensure_output_dir_exists
from plot_choropleth_capacity_factors import plot_choropleth

ROUNDER = 20

MARKET_VALUES = [
    "offwind-ac",
    "offwind-dc",
    "onwind",
    "solar",
    "solar rooftop",
    "ror",
]

MARKET_PRICES = [
    "AC",
    "gas",
    "H2",
    "low voltage",
    "rural heat",
    "urban decentral heat",
    "urban central heat",
]


def get_market_prices(n):
    return (
        n.buses_t.marginal_price.mean()
        .groupby([n.buses.location, n.buses.carrier])
        .first()
        .unstack()
        .drop(["", "EU"])  # drop not spatially resolved
        .dropna(how="all", axis=1)
    )


def get_market_values(n):
    mv_generators = n.statistics.market_value(
        comps={"Generator"}, nice_names=False, groupby=["bus", "carrier"]
    )

    mv_links = n.statistics.market_value(
        comps={"Link"}, nice_names=False, groupby=["bus0", "carrier"]
    ).rename_axis(index={"bus0": "bus"})

    mv = pd.concat([mv_generators, mv_links]).droplevel(0)

    mv.index = pd.MultiIndex.from_arrays(
        [
            mv.index.get_level_values("bus").map(n.buses.location),
            mv.index.get_level_values("carrier"),
        ]
    )

    to_drop = list(n.buses.index[n.buses.index.str.len() == 2]) + ["", "EU", "process"]
    mv = mv.drop(to_drop, errors='ignore').unstack().dropna(how="all", axis=1)

    return mv


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_choropleth_prices",
            simpl="",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    ensure_output_dir_exists(snakemake)

    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
    regions_offshore = gpd.read_file(snakemake.input.regions_offshore).set_index("name")

    n = pypsa.Network(snakemake.input.network)

    lmps = get_market_prices(n)

    vmins = np.floor(lmps.div(ROUNDER)).min().mul(ROUNDER)
    vmaxs = np.ceil(lmps.div(ROUNDER)).max().mul(ROUNDER)

    for carrier in lmps.columns.intersection(MARKET_PRICES):
        plot_choropleth(
            lmps,
            regions_onshore,
            carrier,
            cmap="Spectral_r",
            vmax=vmaxs[carrier],
            vmin=vmins[carrier],
            label="average market price [€/MWh]",
            title=n.carriers.at[carrier, "nice_name"],
            dir=snakemake.output[0] + "/market-prices-",
        )

    mv = get_market_values(n)

    vmins = np.floor(mv.div(ROUNDER)).min().mul(ROUNDER)
    vmaxs = np.ceil(mv.div(ROUNDER)).max().mul(ROUNDER)

    for carrier in mv.columns.intersection(MARKET_VALUES):
        regions = (
            regions_offshore
            if carrier in ["offwind-ac", "offwind-dc"]
            else regions_onshore
        )

        plot_choropleth(
            mv,
            regions,
            carrier,
            cmap="Spectral_r",
            vmax=vmaxs[carrier],
            vmin=vmins[carrier],
            label="average market value [€/MWh]",
            title=n.carriers.at[carrier, "nice_name"],
            dir=snakemake.output[0] + "/market-values-",
        )
