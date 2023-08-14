# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot unclustered electricity transmission network.
"""

import pypsa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.lines import Line2D


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("plot_power_network_topology")

    plt.style.use(snakemake.input.rc)

    n = pypsa.Network(snakemake.input.network)

    w = n.lines.v_nom.div(380)
    c = n.lines.v_nom.map({220: "teal", 300: "orange", 380: "firebrick"})

    fig, ax = plt.subplots(figsize=(13, 13), subplot_kw={"projection": ccrs.EqualEarth()})

    n.plot(
        ax=ax,
        bus_sizes=0,
        line_widths=w,
        line_colors=c,
        link_colors="royalblue",
        link_widths=1,
    )

    handles = [
        Line2D([0], [0], color="teal", lw=2),
        Line2D([0], [0], color="orange", lw=2),
        Line2D([0], [0], color="firebrick", lw=2),
        Line2D([0], [0], color="royalblue", lw=2),
    ]
    plt.legend(
        handles,
        ["HVAC 220 kV", "HVAC 300 kV", "HVAC 380 kV", "HVDC"],
        frameon=False,
        loc=[0.2, 0.85],
        fontsize=14,
    )

    for fn in snakemake.output:
        plt.savefig(fn)
