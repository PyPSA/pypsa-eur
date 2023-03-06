# -*- coding: utf-8 -*-
<<<<<<< HEAD
# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Plots map with pie charts and cost box bar charts.

Relevant Settings
-----------------

Inputs
------

Outputs
-------

Description
-----------
"""

import logging

import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _helpers import (
    aggregate_costs,
    aggregate_p,
    configure_logging,
    load_network_for_plots,
)
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Circle, Ellipse

to_rgba = mpl.colors.colorConverter.to_rgba

logger = logging.getLogger(__name__)


def make_handler_map_to_scale_circles_as_in(ax, dont_resize_actively=False):
    fig = ax.get_figure()

    def axes2pt():
        return np.diff(ax.transData.transform([(0, 0), (1, 1)]), axis=0)[0] * (
            72.0 / fig.dpi
        )

    ellipses = []
    if not dont_resize_actively:

        def update_width_height(event):
            dist = axes2pt()
            for e, radius in ellipses:
                e.width, e.height = 2.0 * radius * dist

        fig.canvas.mpl_connect("resize_event", update_width_height)
        ax.callbacks.connect("xlim_changed", update_width_height)
        ax.callbacks.connect("ylim_changed", update_width_height)

    def legend_circle_handler(
        legend, orig_handle, xdescent, ydescent, width, height, fontsize
    ):
        w, h = 2.0 * orig_handle.get_radius() * axes2pt()
        e = Ellipse(
            xy=(0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent),
            width=w,
            height=w,
        )
        ellipses.append((e, orig_handle.get_radius()))
        return e

    return {Circle: HandlerPatch(patch_func=legend_circle_handler)}


def make_legend_circles_for(sizes, scale=1.0, **kw):
    return [Circle((0, 0), radius=(s / scale) ** 0.5, **kw) for s in sizes]


def set_plot_style():
    plt.style.use(
        [
            "classic",
            "seaborn-white",
            {
                "axes.grid": False,
                "grid.linestyle": "--",
                "grid.color": "0.6",
                "hatch.color": "white",
                "patch.linewidth": 0.5,
                "font.size": 12,
                "legend.fontsize": "medium",
                "lines.linewidth": 1.5,
                "pdf.fonttype": 42,
            },
        ]
    )


def plot_map(n, opts, ax=None, attribute="p_nom"):
    if ax is None:
        ax = plt.gca()

    ## DATA
    line_colors = {
        "cur": "purple",
        "exp": mpl.colors.rgb2hex(to_rgba("red", 0.7), True),
    }
    tech_colors = opts["tech_colors"]

    if attribute == "p_nom":
        # bus_sizes = n.generators_t.p.sum().loc[n.generators.carrier == "load"].groupby(n.generators.bus).sum()
        bus_sizes = pd.concat(
            (
                n.generators.query('carrier != "load"')
                .groupby(["bus", "carrier"])
                .p_nom_opt.sum(),
                n.storage_units.groupby(["bus", "carrier"]).p_nom_opt.sum(),
            )
        )
        line_widths_exp = n.lines.s_nom_opt
        line_widths_cur = n.lines.s_nom_min
        link_widths_exp = n.links.p_nom_opt
        link_widths_cur = n.links.p_nom_min
    else:
        raise "plotting of {} has not been implemented yet".format(attribute)

    line_colors_with_alpha = (line_widths_cur / n.lines.s_nom > 1e-3).map(
        {True: line_colors["cur"], False: to_rgba(line_colors["cur"], 0.0)}
    )
    link_colors_with_alpha = (link_widths_cur / n.links.p_nom > 1e-3).map(
        {True: line_colors["cur"], False: to_rgba(line_colors["cur"], 0.0)}
    )

    ## FORMAT
    linewidth_factor = opts["map"][attribute]["linewidth_factor"]
    bus_size_factor = opts["map"][attribute]["bus_size_factor"]

    ## PLOT
    n.plot(
        line_widths=line_widths_exp / linewidth_factor,
        link_widths=link_widths_exp / linewidth_factor,
        line_colors=line_colors["exp"],
        link_colors=line_colors["exp"],
        bus_sizes=bus_sizes / bus_size_factor,
        bus_colors=tech_colors,
        boundaries=map_boundaries,
        color_geomap=True,
        geomap=True,
        ax=ax,
    )
    n.plot(
        line_widths=line_widths_cur / linewidth_factor,
        link_widths=link_widths_cur / linewidth_factor,
        line_colors=line_colors_with_alpha,
        link_colors=link_colors_with_alpha,
        bus_sizes=0,
        boundaries=map_boundaries,
        color_geomap=True,
        geomap=True,
        ax=ax,
    )
    ax.set_aspect("equal")
    ax.axis("off")

    # Rasterize basemap
    # TODO : Check if this also works with cartopy
    for c in ax.collections[:2]:
        c.set_rasterized(True)

    # LEGEND
    handles = []
    labels = []

    for s in (10, 1):
        handles.append(
            plt.Line2D(
                [0], [0], color=line_colors["exp"], linewidth=s * 1e3 / linewidth_factor
            )
        )
        labels.append("{} GW".format(s))
=======
import logging

logger = logging.getLogger(__name__)

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from helper import override_component_attrs
from make_summary import assign_carriers
from plot_summary import preferred_order, rename_techs
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

plt.style.use(["ggplot", "matplotlibrc"])


def rename_techs_tyndp(tech):
    tech = rename_techs(tech)
    if "heat pump" in tech or "resistive heater" in tech:
        return "power-to-heat"
    elif tech in ["H2 Electrolysis", "methanation", "helmeth", "H2 liquefaction"]:
        return "power-to-gas"
    elif tech == "H2":
        return "H2 storage"
    elif tech in ["NH3", "Haber-Bosch", "ammonia cracker", "ammonia store"]:
        return "ammonia"
    elif tech in ["OCGT", "CHP", "gas boiler", "H2 Fuel Cell"]:
        return "gas-to-power/heat"
    # elif "solar" in tech:
    #     return "solar"
    elif tech in ["Fischer-Tropsch", "methanolisation"]:
        return "power-to-liquid"
    elif "offshore wind" in tech:
        return "offshore wind"
    elif "CC" in tech or "sequestration" in tech:
        return "CCS"
    else:
        return tech


def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = names.str[:i]


def plot_map(
    network,
    components=["links", "stores", "storage_units", "generators"],
    bus_size_factor=1.7e10,
    transmission=False,
    with_legend=True,
):
    tech_colors = snakemake.config["plotting"]["tech_colors"]

    n = network.copy()
    assign_location(n)
    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    costs = pd.DataFrame(index=n.buses.index)

    for comp in components:
        df_c = getattr(n, comp)

        if df_c.empty:
            continue

        df_c["nice_group"] = df_c.carrier.map(rename_techs_tyndp)

        attr = "e_nom_opt" if comp == "stores" else "p_nom_opt"

        costs_c = (
            (df_c.capital_cost * df_c[attr])
            .groupby([df_c.location, df_c.nice_group])
            .sum()
            .unstack()
            .fillna(0.0)
        )
        costs = pd.concat([costs, costs_c], axis=1)

        logger.debug(f"{comp}, {costs}")

    costs = costs.groupby(costs.columns, axis=1).sum()

    costs.drop(list(costs.columns[(costs == 0.0).all()]), axis=1, inplace=True)

    new_columns = preferred_order.intersection(costs.columns).append(
        costs.columns.difference(preferred_order)
    )
    costs = costs[new_columns]

    for item in new_columns:
        if item not in tech_colors:
            logger.warning(f"{item} not in config/plotting/tech_colors")

    costs = costs.stack()  # .sort_index()

    # hack because impossible to drop buses...
    eu_location = snakemake.config["plotting"].get(
        "eu_node_location", dict(x=-5.5, y=46)
    )
    n.buses.loc["EU gas", "x"] = eu_location["x"]
    n.buses.loc["EU gas", "y"] = eu_location["y"]

    n.links.drop(
        n.links.index[(n.links.carrier != "DC") & (n.links.carrier != "B2B")],
        inplace=True,
    )

    # drop non-bus
    to_drop = costs.index.levels[0].symmetric_difference(n.buses.index)
    if len(to_drop) != 0:
        logger.info(f"Dropping non-buses {to_drop.tolist()}")
        costs.drop(to_drop, level=0, inplace=True, axis=0, errors="ignore")

    # make sure they are removed from index
    costs.index = pd.MultiIndex.from_tuples(costs.index.values)

    threshold = 100e6  # 100 mEUR/a
    carriers = costs.groupby(level=1).sum()
    carriers = carriers.where(carriers > threshold).dropna()
    carriers = list(carriers.index)

    # PDF has minimum width, so set these to zero
    line_lower_threshold = 500.0
    line_upper_threshold = 1e4
    linewidth_factor = 4e3
    ac_color = "rosybrown"
    dc_color = "darkseagreen"

    if snakemake.wildcards["lv"] == "1.0":
        # should be zero
        line_widths = n.lines.s_nom_opt - n.lines.s_nom
        link_widths = n.links.p_nom_opt - n.links.p_nom
        title = "added grid"

        if transmission:
            line_widths = n.lines.s_nom_opt
            link_widths = n.links.p_nom_opt
            linewidth_factor = 2e3
            line_lower_threshold = 0.0
            title = "current grid"
    else:
        line_widths = n.lines.s_nom_opt - n.lines.s_nom_min
        link_widths = n.links.p_nom_opt - n.links.p_nom_min
        title = "added grid"

        if transmission:
            line_widths = n.lines.s_nom_opt
            link_widths = n.links.p_nom_opt
            title = "total grid"

    line_widths = line_widths.clip(line_lower_threshold, line_upper_threshold)
    link_widths = link_widths.clip(line_lower_threshold, line_upper_threshold)

    line_widths = line_widths.replace(line_lower_threshold, 0)
    link_widths = link_widths.replace(line_lower_threshold, 0)

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()})
    fig.set_size_inches(7, 6)

    n.plot(
        bus_sizes=costs / bus_size_factor,
        bus_colors=tech_colors,
        line_colors=ac_color,
        link_colors=dc_color,
        line_widths=line_widths / linewidth_factor,
        link_widths=link_widths / linewidth_factor,
        ax=ax,
        **map_opts,
    )

    sizes = [20, 10, 5]
    labels = [f"{s} bEUR/a" for s in sizes]
    sizes = [s / bus_size_factor * 1e9 for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.01, 1.06),
        labelspacing=0.8,
        frameon=False,
        handletextpad=0,
        title="system cost",
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    sizes = [10, 5]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.27, 1.06),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1,
        title=title,
    )

    add_legend_lines(
        ax, sizes, labels, patch_kw=dict(color="lightgrey"), legend_kw=legend_kw
    )

    legend_kw = dict(
        bbox_to_anchor=(1.52, 1.04),
        frameon=False,
    )

    if with_legend:
        colors = [tech_colors[c] for c in carriers] + [ac_color, dc_color]
        labels = carriers + ["HVAC line", "HVDC link"]

        add_legend_patches(
            ax,
            colors,
            labels,
            legend_kw=legend_kw,
        )

    fig.savefig(snakemake.output.map, transparent=True, bbox_inches="tight")


def group_pipes(df, drop_direction=False):
    """
    Group pipes which connect same buses and return overall capacity.
    """
    if drop_direction:
        positive_order = df.bus0 < df.bus1
        df_p = df[positive_order]
        swap_buses = {"bus0": "bus1", "bus1": "bus0"}
        df_n = df[~positive_order].rename(columns=swap_buses)
        df = pd.concat([df_p, df_n])

    # there are pipes for each investment period rename to AC buses name for plotting
    df.index = df.apply(
        lambda x: f"H2 pipeline {x.bus0.replace(' H2', '')} -> {x.bus1.replace(' H2', '')}",
        axis=1,
    )
    # group pipe lines connecting the same buses and rename them for plotting
    pipe_capacity = df.groupby(level=0).agg(
        {"p_nom_opt": sum, "bus0": "first", "bus1": "first"}
    )

    return pipe_capacity


def plot_h2_map(network, regions):
    n = network.copy()
    if "H2 pipeline" not in n.links.carrier.unique():
        return

    assign_location(n)

    h2_storage = n.stores.query("carrier == 'H2'")
    regions["H2"] = h2_storage.rename(
        index=h2_storage.bus.map(n.buses.location)
    ).e_nom_opt.div(
        1e6
    )  # TWh
    regions["H2"] = regions["H2"].where(regions["H2"] > 0.1)

    bus_size_factor = 1e5
    linewidth_factor = 7e3
    # MW below which not drawn
    line_lower_threshold = 750

    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    carriers = ["H2 Electrolysis", "H2 Fuel Cell"]

    elec = n.links[n.links.carrier.isin(carriers)].index

    bus_sizes = (
        n.links.loc[elec, "p_nom_opt"].groupby([n.links["bus0"], n.links.carrier]).sum()
        / bus_size_factor
    )

    # make a fake MultiIndex so that area is correct for legend
    bus_sizes.rename(index=lambda x: x.replace(" H2", ""), level=0, inplace=True)
    # drop all links which are not H2 pipelines
    n.links.drop(
        n.links.index[~n.links.carrier.str.contains("H2 pipeline")], inplace=True
    )

    h2_new = n.links[n.links.carrier == "H2 pipeline"]
    h2_retro = n.links[n.links.carrier == "H2 pipeline retrofitted"]

    if snakemake.config["foresight"] == "myopic":
        # sum capacitiy for pipelines from different investment periods
        h2_new = group_pipes(h2_new)

        if not h2_retro.empty:
            h2_retro = (
                group_pipes(h2_retro, drop_direction=True)
                .reindex(h2_new.index)
                .fillna(0)
            )

    if not h2_retro.empty:
        positive_order = h2_retro.bus0 < h2_retro.bus1
        h2_retro_p = h2_retro[positive_order]
        swap_buses = {"bus0": "bus1", "bus1": "bus0"}
        h2_retro_n = h2_retro[~positive_order].rename(columns=swap_buses)
        h2_retro = pd.concat([h2_retro_p, h2_retro_n])

        h2_retro["index_orig"] = h2_retro.index
        h2_retro.index = h2_retro.apply(
            lambda x: f"H2 pipeline {x.bus0.replace(' H2', '')} -> {x.bus1.replace(' H2', '')}",
            axis=1,
        )

        retro_w_new_i = h2_retro.index.intersection(h2_new.index)
        h2_retro_w_new = h2_retro.loc[retro_w_new_i]

        retro_wo_new_i = h2_retro.index.difference(h2_new.index)
        h2_retro_wo_new = h2_retro.loc[retro_wo_new_i]
        h2_retro_wo_new.index = h2_retro_wo_new.index_orig

        to_concat = [h2_new, h2_retro_w_new, h2_retro_wo_new]
        h2_total = pd.concat(to_concat).p_nom_opt.groupby(level=0).sum()

    else:
        h2_total = h2_new.p_nom_opt

    link_widths_total = h2_total / linewidth_factor

    n.links.rename(index=lambda x: x.split("-2")[0], inplace=True)
    n.links = n.links.groupby(level=0).first()
    link_widths_total = link_widths_total.reindex(n.links.index).fillna(0.0)
    link_widths_total[n.links.p_nom_opt < line_lower_threshold] = 0.0

    retro = n.links.p_nom_opt.where(
        n.links.carrier == "H2 pipeline retrofitted", other=0.0
    )
    link_widths_retro = retro / linewidth_factor
    link_widths_retro[n.links.p_nom_opt < line_lower_threshold] = 0.0

    n.links.bus0 = n.links.bus0.str.replace(" H2", "")
    n.links.bus1 = n.links.bus1.str.replace(" H2", "")

    proj = ccrs.EqualEarth()
    regions = regions.to_crs(proj.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": proj})

    color_h2_pipe = "#b3f3f4"
    color_retrofit = "#499a9c"

    bus_colors = {"H2 Electrolysis": "#ff29d9", "H2 Fuel Cell": "#805394"}

    n.plot(
        geomap=True,
        bus_sizes=bus_sizes,
        bus_colors=bus_colors,
        link_colors=color_h2_pipe,
        link_widths=link_widths_total,
        branch_components=["Link"],
        ax=ax,
        **map_opts,
    )

    n.plot(
        geomap=True,
        bus_sizes=0,
        link_colors=color_retrofit,
        link_widths=link_widths_retro,
        branch_components=["Link"],
        ax=ax,
        color_geomap=False,
        boundaries=map_opts["boundaries"],
    )

    regions.plot(
        ax=ax,
        column="H2",
        cmap="Blues",
        linewidths=0,
        legend=True,
        vmax=6,
        vmin=0,
        legend_kwds={
            "label": "Hydrogen Storage [TWh]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    sizes = [50, 10]
    labels = [f"{s} GW" for s in sizes]
    sizes = [s / bus_size_factor * 1e3 for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1),
        labelspacing=0.8,
        handletextpad=0,
        frameon=False,
    )

    add_legend_circles(
        ax,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    sizes = [30, 10]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.23, 1),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1,
    )

    add_legend_lines(
        ax,
        sizes,
        labels,
        patch_kw=dict(color="lightgrey"),
        legend_kw=legend_kw,
    )

    colors = [bus_colors[c] for c in carriers] + [color_h2_pipe, color_retrofit]
    labels = carriers + ["H2 pipeline (total)", "H2 pipeline (repurposed)"]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0, 1.13),
        ncol=2,
        frameon=False,
    )

    add_legend_patches(ax, colors, labels, legend_kw=legend_kw)

    ax.set_facecolor("white")

    fig.savefig(
        snakemake.output.map.replace("-costs-all", "-h2_network"), bbox_inches="tight"
    )


def plot_ch4_map(network):
    n = network.copy()

    if "gas pipeline" not in n.links.carrier.unique():
        return

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

    methanation_i = n.links[n.links.carrier.isin(["helmeth", "Sabatier"])].index
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

    max_usage = n.links_t.p0.abs().max(axis=0)
    link_widths_used = max_usage / linewidth_factor
    link_widths_used[max_usage < line_lower_threshold] = 0.0

    tech_colors = snakemake.config["plotting"]["tech_colors"]

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

    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": ccrs.EqualEarth()})

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

    fig.savefig(
        snakemake.output.map.replace("-costs-all", "-ch4_network"), bbox_inches="tight"
    )


def plot_map_without(network):
    n = network.copy()
    assign_location(n)

    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": ccrs.EqualEarth()})

    # PDF has minimum width, so set these to zero
    line_lower_threshold = 200.0
    line_upper_threshold = 1e4
    linewidth_factor = 3e3
    ac_color = "rosybrown"
    dc_color = "darkseagreen"

    # hack because impossible to drop buses...
    if "EU gas" in n.buses.index:
        eu_location = snakemake.config["plotting"].get(
            "eu_node_location", dict(x=-5.5, y=46)
        )
        n.buses.loc["EU gas", "x"] = eu_location["x"]
        n.buses.loc["EU gas", "y"] = eu_location["y"]

    to_drop = n.links.index[(n.links.carrier != "DC") & (n.links.carrier != "B2B")]
    n.links.drop(to_drop, inplace=True)

    if snakemake.wildcards["lv"] == "1.0":
        line_widths = n.lines.s_nom
        link_widths = n.links.p_nom
    else:
        line_widths = n.lines.s_nom_min
        link_widths = n.links.p_nom_min

    line_widths = line_widths.clip(line_lower_threshold, line_upper_threshold)
    link_widths = link_widths.clip(line_lower_threshold, line_upper_threshold)

    line_widths = line_widths.replace(line_lower_threshold, 0)
    link_widths = link_widths.replace(line_lower_threshold, 0)

    n.plot(
        bus_colors="k",
        line_colors=ac_color,
        link_colors=dc_color,
        line_widths=line_widths / linewidth_factor,
        link_widths=link_widths / linewidth_factor,
        ax=ax,
        **map_opts,
    )

    handles = []
    labels = []

    for s in (10, 5):
        handles.append(
            plt.Line2D([0], [0], color=ac_color, linewidth=s * 1e3 / linewidth_factor)
        )
        labels.append(f"{s} GW")
>>>>>>> pypsa-eur-sec/master
    l1_1 = ax.legend(
        handles,
        labels,
        loc="upper left",
<<<<<<< HEAD
        bbox_to_anchor=(0.24, 1.01),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1.5,
        title="Transmission Exp./Exist.             ",
    )
    ax.add_artist(l1_1)

    handles = []
    labels = []
    for s in (10, 5):
        handles.append(
            plt.Line2D(
                [0], [0], color=line_colors["cur"], linewidth=s * 1e3 / linewidth_factor
            )
        )
        labels.append("/")
    l1_2 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.26, 1.01),
        frameon=False,
        labelspacing=0.8,
        handletextpad=0.5,
        title=" ",
    )
    ax.add_artist(l1_2)

    handles = make_legend_circles_for(
        [10e3, 5e3, 1e3], scale=bus_size_factor, facecolor="w"
    )
    labels = ["{} GW".format(s) for s in (10, 5, 3)]
    l2 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.01, 1.01),
        frameon=False,
        labelspacing=1.0,
        title="Generation",
        handler_map=make_handler_map_to_scale_circles_as_in(ax),
    )
    ax.add_artist(l2)

    techs = (bus_sizes.index.levels[1]).intersection(
        pd.Index(opts["vre_techs"] + opts["conv_techs"] + opts["storage_techs"])
    )
    handles = []
    labels = []
    for t in techs:
        handles.append(
            plt.Line2D(
                [0], [0], color=tech_colors[t], marker="o", markersize=8, linewidth=0
            )
        )
        labels.append(opts["nice_names"].get(t, t))
    l3 = ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.0),  # bbox_to_anchor=(0.72, -0.05),
        handletextpad=0.0,
        columnspacing=0.5,
        ncol=4,
        title="Technology",
    )

    return fig


def plot_total_energy_pie(n, opts, ax=None):
    if ax is None:
        ax = plt.gca()

    ax.set_title("Energy per technology", fontdict=dict(fontsize="medium"))

    e_primary = aggregate_p(n).drop("load", errors="ignore").loc[lambda s: s > 0]

    patches, texts, autotexts = ax.pie(
        e_primary,
        startangle=90,
        labels=e_primary.rename(opts["nice_names"]).index,
        autopct="%.0f%%",
        shadow=False,
        colors=[opts["tech_colors"][tech] for tech in e_primary.index],
    )
    for t1, t2, i in zip(texts, autotexts, e_primary.index):
        if e_primary.at[i] < 0.04 * e_primary.sum():
            t1.remove()
            t2.remove()


def plot_total_cost_bar(n, opts, ax=None):
    if ax is None:
        ax = plt.gca()

    total_load = (n.snapshot_weightings.generators * n.loads_t.p.sum(axis=1)).sum()
    tech_colors = opts["tech_colors"]

    def split_costs(n):
        costs = aggregate_costs(n).reset_index(level=0, drop=True)
        costs_ex = aggregate_costs(n, existing_only=True).reset_index(
            level=0, drop=True
        )
        return (
            costs["capital"].add(costs["marginal"], fill_value=0.0),
            costs_ex["capital"],
            costs["capital"] - costs_ex["capital"],
            costs["marginal"],
        )

    costs, costs_cap_ex, costs_cap_new, costs_marg = split_costs(n)

    costs_graph = pd.DataFrame(
        dict(a=costs.drop("load", errors="ignore")),
        index=[
            "AC-AC",
            "AC line",
            "onwind",
            "offwind-ac",
            "offwind-dc",
            "solar",
            "OCGT",
            "CCGT",
            "battery",
            "H2",
        ],
    ).dropna()
    bottom = np.array([0.0, 0.0])
    texts = []

    for i, ind in enumerate(costs_graph.index):
        data = np.asarray(costs_graph.loc[ind]) / total_load
        ax.bar([0.5], data, bottom=bottom, color=tech_colors[ind], width=0.7, zorder=-1)
        bottom_sub = bottom
        bottom = bottom + data

        if ind in opts["conv_techs"] + ["AC line"]:
            for c in [costs_cap_ex, costs_marg]:
                if ind in c:
                    data_sub = np.asarray([c.loc[ind]]) / total_load
                    ax.bar(
                        [0.5],
                        data_sub,
                        linewidth=0,
                        bottom=bottom_sub,
                        color=tech_colors[ind],
                        width=0.7,
                        zorder=-1,
                        alpha=0.8,
                    )
                    bottom_sub += data_sub

        if abs(data[-1]) < 5:
            continue

        text = ax.text(
            1.1, (bottom - 0.5 * data)[-1] - 3, opts["nice_names"].get(ind, ind)
        )
        texts.append(text)

    ax.set_ylabel("Average system cost [Eur/MWh]")
    ax.set_ylim([0, opts.get("costs_max", 80)])
    ax.set_xlim([0, 1])
    ax.set_xticklabels([])
    ax.grid(True, axis="y", color="k", linestyle="dotted")
=======
        bbox_to_anchor=(0.05, 1.01),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1.5,
        title="Today's transmission",
    )
    ax.add_artist(l1_1)

    fig.savefig(snakemake.output.today, transparent=True, bbox_inches="tight")


def plot_series(network, carrier="AC", name="test"):
    n = network.copy()
    assign_location(n)
    assign_carriers(n)

    buses = n.buses.index[n.buses.carrier.str.contains(carrier)]

    supply = pd.DataFrame(index=n.snapshots)
    for c in n.iterate_components(n.branch_components):
        n_port = 4 if c.name == "Link" else 2
        for i in range(n_port):
            supply = pd.concat(
                (
                    supply,
                    (-1)
                    * c.pnl["p" + str(i)]
                    .loc[:, c.df.index[c.df["bus" + str(i)].isin(buses)]]
                    .groupby(c.df.carrier, axis=1)
                    .sum(),
                ),
                axis=1,
            )

    for c in n.iterate_components(n.one_port_components):
        comps = c.df.index[c.df.bus.isin(buses)]
        supply = pd.concat(
            (
                supply,
                ((c.pnl["p"].loc[:, comps]).multiply(c.df.loc[comps, "sign"]))
                .groupby(c.df.carrier, axis=1)
                .sum(),
            ),
            axis=1,
        )

    supply = supply.groupby(rename_techs_tyndp, axis=1).sum()

    both = supply.columns[(supply < 0.0).any() & (supply > 0.0).any()]

    positive_supply = supply[both]
    negative_supply = supply[both]

    positive_supply[positive_supply < 0.0] = 0.0
    negative_supply[negative_supply > 0.0] = 0.0

    supply[both] = positive_supply

    suffix = " charging"

    negative_supply.columns = negative_supply.columns + suffix

    supply = pd.concat((supply, negative_supply), axis=1)

    # 14-21.2 for flaute
    # 19-26.1 for flaute

    start = "2013-02-19"
    stop = "2013-02-26"

    threshold = 10e3

    to_drop = supply.columns[(abs(supply) < threshold).all()]

    if len(to_drop) != 0:
        logger.info(f"Dropping {to_drop.tolist()} from supply")
        supply.drop(columns=to_drop, inplace=True)

    supply.index.name = None

    supply = supply / 1e3

    supply.rename(
        columns={"electricity": "electric demand", "heat": "heat demand"}, inplace=True
    )
    supply.columns = supply.columns.str.replace("residential ", "")
    supply.columns = supply.columns.str.replace("services ", "")
    supply.columns = supply.columns.str.replace("urban decentral ", "decentral ")

    preferred_order = pd.Index(
        [
            "electric demand",
            "transmission lines",
            "hydroelectricity",
            "hydro reservoir",
            "run of river",
            "pumped hydro storage",
            "CHP",
            "onshore wind",
            "offshore wind",
            "solar PV",
            "solar thermal",
            "building retrofitting",
            "ground heat pump",
            "air heat pump",
            "resistive heater",
            "OCGT",
            "gas boiler",
            "gas",
            "natural gas",
            "methanation",
            "hydrogen storage",
            "battery storage",
            "hot water storage",
        ]
    )

    new_columns = preferred_order.intersection(supply.columns).append(
        supply.columns.difference(preferred_order)
    )

    supply = supply.groupby(supply.columns, axis=1).sum()
    fig, ax = plt.subplots()
    fig.set_size_inches((8, 5))

    (
        supply.loc[start:stop, new_columns].plot(
            ax=ax,
            kind="area",
            stacked=True,
            linewidth=0.0,
            color=[
                snakemake.config["plotting"]["tech_colors"][i.replace(suffix, "")]
                for i in new_columns
            ],
        )
    )

    handles, labels = ax.get_legend_handles_labels()

    handles.reverse()
    labels.reverse()

    new_handles = []
    new_labels = []

    for i, item in enumerate(labels):
        if "charging" not in item:
            new_handles.append(handles[i])
            new_labels.append(labels[i])

    ax.legend(new_handles, new_labels, ncol=3, loc="upper left", frameon=False)
    ax.set_xlim([start, stop])
    ax.set_ylim([-1300, 1900])
    ax.grid(True)
    ax.set_ylabel("Power [GW]")
    fig.tight_layout()

    fig.savefig(
        "{}{}/maps/series-{}-{}-{}-{}-{}.pdf".format(
            snakemake.config["results_dir"],
            snakemake.config["run"],
            snakemake.wildcards["lv"],
            carrier,
            start,
            stop,
            name,
        ),
        transparent=True,
    )
>>>>>>> pypsa-eur-sec/master


if __name__ == "__main__":
    if "snakemake" not in globals():
<<<<<<< HEAD
        from _helpers import mock_snakemake
=======
        from helper import mock_snakemake
>>>>>>> pypsa-eur-sec/master

        snakemake = mock_snakemake(
            "plot_network",
            simpl="",
<<<<<<< HEAD
            clusters="5",
            ll="copt",
            opts="Co2L-24H",
            attr="p_nom",
            ext="pdf",
        )
    configure_logging(snakemake)

    set_plot_style()

    config, wildcards = snakemake.config, snakemake.wildcards

    map_figsize = config["plotting"]["map"]["figsize"]
    map_boundaries = config["plotting"]["map"]["boundaries"]

    n = load_network_for_plots(
        snakemake.input.network, snakemake.input.tech_costs, config
    )

    scenario_opts = wildcards.opts.split("-")

    fig, ax = plt.subplots(
        figsize=map_figsize, subplot_kw={"projection": ccrs.PlateCarree()}
    )
    plot_map(n, config["plotting"], ax=ax, attribute=wildcards.attr)

    fig.savefig(snakemake.output.only_map, dpi=150, bbox_inches="tight")

    ax1 = fig.add_axes([-0.115, 0.625, 0.2, 0.2])
    plot_total_energy_pie(n, config["plotting"], ax=ax1)

    ax2 = fig.add_axes([-0.075, 0.1, 0.1, 0.45])
    plot_total_cost_bar(n, config["plotting"], ax=ax2)

    ll = wildcards.ll
    ll_type = ll[0]
    ll_factor = ll[1:]
    lbl = dict(c="line cost", v="line volume")[ll_type]
    amnt = "{ll} x today's".format(ll=ll_factor) if ll_factor != "opt" else "optimal"
    fig.suptitle(
        "Expansion to {amount} {label} at {clusters} clusters".format(
            amount=amnt, label=lbl, clusters=wildcards.clusters
        )
    )

    fig.savefig(snakemake.output.ext, transparent=True, bbox_inches="tight")
=======
            clusters="181",
            lv="opt",
            opts="",
            sector_opts="Co2L0-730H-T-H-B-I-A-solar+p3-linemaxext10",
            planning_horizons="2050",
        )

    logging.basicConfig(level=snakemake.config["logging_level"])

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    regions = gpd.read_file(snakemake.input.regions).set_index("name")

    map_opts = snakemake.config["plotting"]["map"]

    plot_map(
        n,
        components=["generators", "links", "stores", "storage_units"],
        bus_size_factor=2e10,
        transmission=False,
    )

    plot_h2_map(n, regions)
    plot_ch4_map(n)
    plot_map_without(n)

    # plot_series(n, carrier="AC", name=suffix)
    # plot_series(n, carrier="heat", name=suffix)
>>>>>>> pypsa-eur-sec/master
