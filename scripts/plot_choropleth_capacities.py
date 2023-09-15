# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot optimised capacities on map.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from plot_choropleth_capacity_factors import plot_choropleth

ROUNDER = 10

IGNORE_LINKS = [
    "DC",
    "H2 pipeline",
    "H2 pipeline retrofitted",
    "gas pipeline",
    "gas pipeline new",
]

CAPACITIES = [
    "Fischer-Tropsch",
    "H2 Electrolysis",
    "H2 Fuel Cell",
    "Haber-Bosch",
    "OCGT",
    "PHS",
    "SMR",
    "SMR CC",
    "Sabatier",
    "battery charger",
    "battery discharger",
    "electricity distribution grid",
    "home battery charger",
    "home battery discharger",
    "hydro",
    "methanolisation",
    "offwind-ac",
    "offwind-dc",
    "onwind",
    "ror",
    "gas boiler",
    "air heat pump",
    "ground heat pump",
    "resistive heater",
    "solar",
    "solar rooftop",
    "gas CHP",
    "gas CHP CC",
]

PREFIXES_TO_REMOVE = [
    "residential ",
    "services ",
    "urban ",
    "rural ",
    "central ",
    "decentral ",
    "home ",
]


def remove_prefixes(s):
    for prefix in PREFIXES_TO_REMOVE:
        s = s.replace(prefix, "", 1)
    return s


def get_optimal_capacity(n):
    p_nom_opt_oneport = n.statistics.optimal_capacity(
        comps={"Generator", "StorageUnit"}, groupby=["bus", "carrier"]
    )
    p_nom_opt_links = n.statistics.optimal_capacity(
        comps={"Link"}, groupby=["bus0", "carrier"]
    ).rename_axis(index={"bus0": "bus"})
    p_nom_opt = (
        pd.concat([p_nom_opt_oneport, p_nom_opt_links]).droplevel(0).div(1e3)
    )  # GW

    p_nom_opt.index = pd.MultiIndex.from_arrays(
        [
            p_nom_opt.index.get_level_values("bus").map(n.buses.location),
            p_nom_opt.index.get_level_values("carrier"),
        ]
    )

    p_nom_opt = p_nom_opt.unstack().drop(["", "EU"]).dropna(how="all", axis=1)

    p_nom_opt = p_nom_opt.loc[:, ~p_nom_opt.columns.isin(IGNORE_LINKS)]

    return p_nom_opt.groupby(p_nom_opt.columns.map(remove_prefixes), axis=1).sum()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_choropleth_capacities",
            simpl="",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
    regions_offshore = gpd.read_file(snakemake.input.regions_offshore).set_index("name")

    n = pypsa.Network(snakemake.input.network)

    p_nom_opt = get_optimal_capacity(n)

    vmins = np.floor(p_nom_opt.div(ROUNDER)).min().mul(ROUNDER)
    vmaxs = np.ceil(p_nom_opt.div(ROUNDER)).max().mul(ROUNDER)

    for carrier in p_nom_opt.columns.intersection(CAPACITIES):
        regions = (
            regions_offshore
            if carrier in ["offwind-ac", "offwind-dc"]
            else regions_onshore
        )

        plot_choropleth(
            p_nom_opt,
            regions,
            carrier,
            cmap="YlGnBu",
            vmax=vmaxs[carrier],
            vmin=vmins[carrier],
            label=f"optimised capacity [GW]",
            title=n.carriers.nice_name.get(carrier, carrier),
            dir=snakemake.output[0],
        )
