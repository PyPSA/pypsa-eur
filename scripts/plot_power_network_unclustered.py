# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot unclustered electricity transmission network.
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pypsa
from matplotlib.lines import Line2D
import seaborn as sns

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_power_network_unclustered",
            configfiles="config/config.20240826-z1.yaml",
        )

    plt.style.use(snakemake.input.rc)

    n = pypsa.Network(snakemake.input.network)

    voltages = n.lines.v_nom.unique()
    cmap = sns.color_palette("rainbow", as_cmap=True)
    norm = plt.Normalize(voltages.min(), voltages.max())
    color_mapping_dict = {v: cmap(norm(v)) for v in voltages}

    w = n.lines.v_nom.div(380 / 1.5)
    c = n.lines.v_nom.map(color_mapping_dict)

    fig, ax = plt.subplots(
        figsize=(13, 13), subplot_kw={"projection": ccrs.EqualEarth()}
    )

    n.plot(
        ax=ax,
        bus_sizes=0,
        line_widths=w,
        line_colors=c,
        link_colors="royalblue",
        link_widths=2,
    )

    handles = [
        Line2D([0], [0], color=color_mapping_dict[220], lw=2),
        Line2D([0], [0], color=color_mapping_dict[275], lw=2),
        Line2D([0], [0], color=color_mapping_dict[300], lw=2),
        Line2D([0], [0], color=color_mapping_dict[330], lw=2),
        Line2D([0], [0], color=color_mapping_dict[380], lw=2),
        Line2D([0], [0], color=color_mapping_dict[400], lw=2),
        Line2D([0], [0], color="royalblue", lw=2),  # HVDC remains royalblue
    ]

    plt.legend(
        handles,
        [
            "HVAC 220 kV",
            "HVAC 275 kV",
            "HVAC 300 kV",
            "HVAC 330 kV",
            "HVAC 380 kV",
            "HVAC 400 kV",
            "HVDC",
        ],
        frameon=False,
        loc=[0.2, 0.85],
        fontsize=14,
    )

    for fn in snakemake.output:
        plt.savefig(fn)
