# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates map of optimised gas network, storage and selected other
infrastructure.
"""

import logging

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from _helpers import configure_logging
from plot_power_network import assign_location, load_projection
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

logger = logging.getLogger(__name__)


def plot_ch4_map(n):
    # if "gas pipeline" not in n.links.carrier.unique():
    #     return

    assign_location(n)

    bus_size_factor = 8e7
    linewidth_factor = 1e4
    # MW below which not drawn
    line_lower_threshold = 1e3

    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    fossil_gas_i = n.generators[n.generators.carrier == "gas"].index
    fossil_gas = (
        n.generators_t.p.loc[:, fossil_gas_i]
        .mul(n.snapshot_weightings.generators, axis=0)
        .sum()
        .groupby(n.generators.loc[fossil_gas_i, "bus"])
        .sum()
        / bus_size_factor
    )
    fossil_gas.rename(index=lambda x: x.replace(" gas", ""), inplace=True)
    fossil_gas = fossil_gas.reindex(n.buses.index).fillna(0)
    # make a fake MultiIndex so that area is correct for legend
    fossil_gas.index = pd.MultiIndex.from_product([fossil_gas.index, ["fossil gas"]])

    methanation_i = n.links.query("carrier == 'Sabatier'").index
    methanation = (
        abs(
            n.links_t.p1.loc[:, methanation_i].mul(
                n.snapshot_weightings.generators, axis=0
            )
        )
        .sum()
        .groupby(n.links.loc[methanation_i, "bus1"])
        .sum()
        / bus_size_factor
    )
    methanation = (
        methanation.groupby(methanation.index)
        .sum()
        .rename(index=lambda x: x.replace(" gas", ""))
    )
    # make a fake MultiIndex so that area is correct for legend
    methanation.index = pd.MultiIndex.from_product([methanation.index, ["methanation"]])

    biogas_i = n.stores[n.stores.carrier == "biogas"].index
    biogas = (
        n.stores_t.p.loc[:, biogas_i]
        .mul(n.snapshot_weightings.generators, axis=0)
        .sum()
        .groupby(n.stores.loc[biogas_i, "bus"])
        .sum()
        / bus_size_factor
    )
    biogas = (
        biogas.groupby(biogas.index)
        .sum()
        .rename(index=lambda x: x.replace(" biogas", ""))
    )
    # make a fake MultiIndex so that area is correct for legend
    biogas.index = pd.MultiIndex.from_product([biogas.index, ["biogas"]])

    bus_sizes = pd.concat([fossil_gas, methanation, biogas])
    bus_sizes.sort_index(inplace=True)

    to_remove = n.links.index[~n.links.carrier.str.contains("gas pipeline")]
    n.links.drop(to_remove, inplace=True)

    link_widths_rem = n.links.p_nom_opt / linewidth_factor
    link_widths_rem[n.links.p_nom_opt < line_lower_threshold] = 0.0

    link_widths_orig = n.links.p_nom / linewidth_factor
    link_widths_orig[n.links.p_nom < line_lower_threshold] = 0.0

    max_usage = n.links_t.p0[n.links.index].abs().max(axis=0)
    link_widths_used = max_usage / linewidth_factor
    link_widths_used[max_usage < line_lower_threshold] = 0.0

    tech_colors = snakemake.params.plotting["tech_colors"]

    pipe_colors = {
        "gas pipeline": "#f08080",
        "gas pipeline new": "#c46868",
        "gas pipeline (in 2020)": "lightgrey",
        "gas pipeline (available)": "#e8d1d1",
    }

    link_color_used = n.links.carrier.map(pipe_colors)

    n.links.bus0 = n.links.bus0.str.replace(" gas", "")
    n.links.bus1 = n.links.bus1.str.replace(" gas", "")

    bus_colors = {
        "fossil gas": tech_colors["fossil gas"],
        "methanation": tech_colors["methanation"],
        "biogas": "seagreen",
    }

    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": proj})

    n.plot(
        bus_sizes=bus_sizes,
        bus_colors=bus_colors,
        link_colors=pipe_colors["gas pipeline (in 2020)"],
        link_widths=link_widths_orig,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    n.plot(
        ax=ax,
        bus_sizes=0.0,
        link_colors=pipe_colors["gas pipeline (available)"],
        link_widths=link_widths_rem,
        branch_components=["Link"],
        color_geomap=False,
        boundaries=map_opts["boundaries"],
    )

    n.plot(
        ax=ax,
        bus_sizes=0.0,
        link_colors=link_color_used,
        link_widths=link_widths_used,
        branch_components=["Link"],
        color_geomap=False,
        boundaries=map_opts["boundaries"],
    )

    sizes = [100, 10]
    labels = [f"{s} TWh" for s in sizes]
    sizes = [s / bus_size_factor * 1e6 for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1.03),
        labelspacing=0.8,
        frameon=False,
        handletextpad=1,
        title="gas sources",
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    sizes = [50, 10]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.25, 1.03),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1,
        title="gas pipeline",
    )

    add_legend_lines(
        ax,
        sizes,
        labels,
        patch_kw=dict(color="lightgrey"),
        legend_kw=legend_kw,
    )

    colors = list(pipe_colors.values()) + list(bus_colors.values())
    labels = list(pipe_colors.keys()) + list(bus_colors.keys())

    # legend on the side
    # legend_kw = dict(
    #     bbox_to_anchor=(1.47, 1.04),
    #     frameon=False,
    # )

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1.24),
        ncol=2,
        frameon=False,
    )

    add_legend_patches(
        ax,
        colors,
        labels,
        legend_kw=legend_kw,
    )

    fig.savefig(snakemake.output.map, bbox_inches="tight")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_gas_network",
            simpl="",
            opts="",
            clusters="37",
            ll="v1.0",
            sector_opts="4380H-T-H-B-I-A-dist1",
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    regions = gpd.read_file(snakemake.input.regions).set_index("name")

    map_opts = snakemake.params.plotting["map"]

    if map_opts["boundaries"] is None:
        map_opts["boundaries"] = regions.total_bounds[[0, 2, 1, 3]] + [-1, 1, -1, 1]

    proj = load_projection(snakemake.params.plotting)

    plot_ch4_map(n)
