# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Create energy balance maps for the defined carriers.
"""
# %%
import logging

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
import seaborn as sns
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches
from pypsa.statistics import get_transmission_carriers

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

    snakemake = mock_snakemake(
        "plot_balance_map",
        simpl="",
        opts="",
        clusters="70",
        ll="vopt",
        sector_opts="",
        planning_horizons="2050",
        run="maps",
        carrier="gas",
    )

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)
    n = pypsa.Network(snakemake.input.network)
    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    plotting = snakemake.params.plotting
    carrier = snakemake.wildcards.carrier
    if carrier not in n.buses.carrier.unique():
        bus_carrier = [c for c in n.buses.carrier.unique() if carrier in c]
        if bus_carrier:
            logging.warning(
                f"Carrier {carrier} is not in the network. Using {bus_carrier} instead."
            )
        else:
            raise ValueError(f"Carrier {carrier} is not in the network.")
    else:
        bus_carrier = carrier

    sns.set_theme(**plotting.get("theme", {}))

    plotting = plotting.get("balance_map", {})
    fig_size = plotting.get("fig_size", (6, 6))
    alpha = plotting.get("alpha", 1)
    region_alpha = plotting.get("region_alpha", 0.6)
    fig, ax = plt.subplots(
        figsize=fig_size,
        subplot_kw={"projection": ccrs.EqualEarth()},
        layout="constrained",
    )

    carrier_plotting = plotting.get(carrier, {})

    n.buses.x = n.buses.location.map(n.buses.x)
    n.buses.y = n.buses.location.map(n.buses.y)

    s = n.statistics
    s.set_parameters(round=3)
    grouper = s.groupers.get_bus_and_carrier
    # %%
    # bus_sizes according to energy balance of bus carrier
    conversion = float(carrier_plotting.get("unit_conversion", 1e6))
    energy_balance_df = s.energy_balance(
        nice_names=True, bus_carrier=bus_carrier, groupby=grouper
    ).round(2)
    # remove energy balance of transmission carriers, which are relate to losses
    transmission_carriers = get_transmission_carriers(
        n, bus_carrier=bus_carrier
    ).rename(
        {"name": "carrier"}
    )  # TODO change in pypsa
    # %%
    energy_balance_df.loc[transmission_carriers.unique(0)] = energy_balance_df.loc[
        transmission_carriers.unique(0)
    ].drop(index=transmission_carriers.unique(1), level="carrier")
    energy_balance_df = energy_balance_df.dropna()
    # %%
    bus_sizes = (
        energy_balance_df.groupby(level=["bus", "carrier"]).sum().div(conversion)
    )

    colors = (
        bus_sizes.index.get_level_values("carrier")
        .unique()
        .to_series()
        .map(n.carriers.set_index("nice_name").color)
    )

    # %%
    # line and links widths according to optimal capacity
    flow = (
        s.transmission(groupby=False, bus_carrier=bus_carrier).div(conversion).round(2)
    )

    if not flow.index.get_level_values(1).empty:
        flow_reversed_mask = flow.index.get_level_values(1).str.contains("reversed")
        flow_reversed = flow[flow_reversed_mask].rename(
            lambda x: x.replace("-reversed", "")
        )
        flow = flow[~flow_reversed_mask].subtract(flow_reversed, fill_value=0)

    # if there are not lines or links for the bus carrier, use fallback for plotting
    fallback = pd.Series()
    line_widths = flow.get("Line", fallback).abs()
    link_widths = flow.get("Link", fallback).abs()

    # %%
    # define maximal size of buses and branch width
    bus_size_factor = float(carrier_plotting.get("bus_factor", 2e-5))
    branch_width_factor = float(carrier_plotting.get("branch_factor", 2e-4))
    flow_size_factor = float(carrier_plotting.get("flow_factor", 2e-4))

    # get prices per region as colormap
    buses = n.buses.query("carrier in @bus_carrier").index
    price = (
        n.buses_t.marginal_price.mean()
        .reindex(buses)
        .rename(n.buses.location)
        .groupby(level=0)
        .mean()
    )
    regions["price"] = price.reindex(regions.index).fillna(0)

    vmin, vmax = regions.price.min(), regions.price.max()
    cmap = carrier_plotting.get("region_cmap", "Greens")

    regions.plot(
        ax=ax,
        column="price",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        edgecolor="None",
        linewidth=0,
        alpha=region_alpha,
        transform=ccrs.PlateCarree(),
        aspect="equal",
    )

    # %%
    # plot map
    n.plot(
        bus_sizes=bus_sizes * bus_size_factor,
        bus_colors=colors,
        bus_alpha=alpha,
        bus_split_circles=True,
        line_widths=link_widths * branch_width_factor,
        link_widths=link_widths * branch_width_factor,
        flow=flow * flow_size_factor,
        ax=ax,
        margin=0.2,
        color_geomap={"border": "darkgrey", "coastline": "darkgrey"},
        geomap="10m",
    )

    # TODO with config
    ax.set_title("Balance Map of carrier " + bus_carrier)

    # Add legend
    # %%
    legend_kwargs = {
        "loc": "upper left",
        "frameon": True,
        "framealpha": 0.5,
        "edgecolor": "None",
    }

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    price_unit = carrier_plotting.get("region_unit", "€/MWh")
    carrier = n.carriers.loc[bus_carrier, "nice_name"]
    cbr = fig.colorbar(
        sm,
        ax=ax,
        label=f"Average Marginal Price {carrier} [{price_unit}]",
        shrink=0.95,
        pad=0.03,
        aspect=50,
        alpha=region_alpha,
        orientation="horizontal",
    )
    cbr.outline.set_edgecolor("None")

    pad = 0.18
    carriers = n.carriers.set_index("nice_name")
    carriers.loc["", "color"] = "None"
    prod_carriers = bus_sizes[bus_sizes > 0].index.unique("carrier").sort_values()
    cons_carriers = (
        bus_sizes[bus_sizes < 0]
        .index.unique("carrier")
        .difference(prod_carriers)
        .sort_values()
    )

    # Add production carriers
    add_legend_patches(
        ax,
        carriers.color[prod_carriers],
        prod_carriers,
        patch_kw={"alpha": alpha},
        legend_kw={
            "bbox_to_anchor": (0, -pad),
            "ncol": 1,
            "title": "Production",
            **legend_kwargs,
        },
    )

    # Add consumption carriers
    add_legend_patches(
        ax,
        carriers.color[cons_carriers],
        cons_carriers,
        patch_kw={"alpha": alpha},
        legend_kw={
            "bbox_to_anchor": (0.5, -pad),
            "ncol": 1,
            "title": "Consumption",
            **legend_kwargs,
        },
    )

    # TODO with config
    # Add bus legend
    legend_bus_sizes = np.array(carrier_plotting.get("bus_sizes", [10, 50]))
    carrier_unit = carrier_plotting.get("unit", "TWh")
    if legend_bus_sizes is not None:
        add_legend_circles(
            ax,
            [s * bus_size_factor for s in legend_bus_sizes],
            [f"{s} {carrier_unit}" for s in legend_bus_sizes],
            legend_kw={
                "bbox_to_anchor": (0, 1),
                "title": "Supply/Demand",
                **legend_kwargs,
            },
        )

    legend_branch_sizes = carrier_plotting.get("branch_sizes", [1, 10])
    if legend_branch_sizes:
        # Add branch legend
        if legend_branch_sizes is not None:
            add_legend_lines(
                ax,
                [s * branch_width_factor for s in legend_branch_sizes],
                [f"{s} {carrier_unit}" for s in legend_branch_sizes],
                legend_kw={"bbox_to_anchor": (0, 0.85), **legend_kwargs},
            )

    fig.savefig(
        snakemake.output.map,
        dpi=300,
        bbox_inches="tight",
    )
