# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot share of potential used on map.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
from _helpers import ensure_output_dir_exists
from plot_choropleth_capacity_factors import plot_choropleth

POTENTIAL = [
    "offwind-ac",
    "offwind-dc",
    "onwind",
    "solar",
    "solar rooftop",
]


def get_potential_used(n):
    return (
        n.generators.eval("p_nom_opt/p_nom_max*100")
        .groupby([n.generators.bus.map(n.buses.location), n.generators.carrier])
        .sum()
        .unstack()
        .drop("EU")
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_choropleth_potential_used",
            simpl="",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    ensure_output_dir_exists(snakemake)

    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
    regions_offshore = gpd.read_file(snakemake.input.regions_offshore).set_index("name")

    n = pypsa.Network(snakemake.input.network)

    potential_used = get_potential_used(n)

    legend_kwds = {
        "label": "potential used [%]",
        "shrink": 0.7,
    }

    for carrier in potential_used.columns.intersection(POTENTIAL):
        regions = (
            regions_offshore
            if carrier in ["offwind-ac", "offwind-dc"]
            else regions_onshore
        )

        cmap = "Reds" if "solar" in carrier else "Blues"

        plot_choropleth(
            potential_used,
            regions,
            carrier,
            cmap=cmap,
            vmax=100,
            vmin=0,
            legend_kwds=legend_kwds,
            title=n.carriers.nice_name.get(carrier, carrier),
            dir=snakemake.output[0],
        )
