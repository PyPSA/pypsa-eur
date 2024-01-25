# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023-2024 PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Plot clustered electricity transmission network.
"""

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
from matplotlib.lines import Line2D
from plot_power_network import load_projection
from pypsa.plot import add_legend_lines

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_power_network_clustered",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )

    lw_factor = 2e3

    n = pypsa.Network(snakemake.input.network)

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    proj = load_projection(snakemake.params.plotting)

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"projection": proj})
    regions.to_crs(proj.proj4_init).plot(
        ax=ax, facecolor="none", edgecolor="lightgray", linewidth=0.75
    )
    n.plot(
        ax=ax,
        margin=0.06,
        line_widths=n.lines.s_nom / lw_factor,
        link_colors=n.links.p_nom.apply(
            lambda x: "darkseagreen" if x > 0 else "skyblue"
        ),
        link_widths=2.0,
    )

    sizes = [10, 20]
    labels = [f"HVAC ({s} GW)" for s in sizes]
    scale = 1e3 / lw_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc=[0.25, 0.9],
        frameon=False,
        labelspacing=0.5,
        handletextpad=1,
        fontsize=13,
    )

    add_legend_lines(
        ax, sizes, labels, patch_kw=dict(color="rosybrown"), legend_kw=legend_kw
    )

    handles = [
        Line2D([0], [0], color="darkseagreen", lw=2),
        Line2D([0], [0], color="skyblue", lw=2),
    ]
    plt.legend(
        handles,
        ["HVDC existing", "HVDC planned"],
        frameon=False,
        loc=[0.0, 0.9],
        fontsize=13,
    )

    plt.savefig(snakemake.output.map, bbox_inches="tight")
