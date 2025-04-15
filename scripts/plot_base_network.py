# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Plot base network transmission network.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
from _helpers import set_scenario_config
from plot_power_network import load_projection
from pypsa.plot import add_legend_lines

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("plot_base_network", run="tyndp-raw")
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.network)

    lw_factor = 1e3 if n.lines.empty else 2e3

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
        link_widths=n.links.p_nom / lw_factor,
    )

    if not n.lines.empty:
        sizes_ac = [10, 20]
        labels_ac = [f"HVAC ({s} GW)" for s in sizes_ac]
        scale_ac = 1e3 / lw_factor
        sizes_ac = [s * scale_ac for s in sizes_ac]

        legend_kw_ac = dict(
            loc=[0.25, 0.9],
            frameon=False,
            labelspacing=0.5,
            handletextpad=1,
            fontsize=13,
        )

        add_legend_lines(
            ax,
            sizes_ac,
            labels_ac,
            patch_kw=dict(color="rosybrown"),
            legend_kw=legend_kw_ac,
        )

    if not n.links.empty:
        sizes_dc = [1, 5]
        labels_dc = [f"HVDC ({s} GW)" for s in sizes_dc]
        scale_dc = 1e3 / lw_factor
        sizes_dc = [s * scale_dc for s in sizes_dc]

        legend_kw_dc = dict(
            loc=[0.0, 0.9],
            frameon=False,
            labelspacing=0.5,
            handletextpad=1,
            fontsize=13,
        )

        add_legend_lines(
            ax,
            sizes_dc,
            labels_dc,
            patch_kw=dict(color="darkseagreen"),
            legend_kw=legend_kw_dc,
        )

    plt.savefig(snakemake.output.map, bbox_inches="tight")
    plt.close()
