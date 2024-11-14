# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023-2024 PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Plot map of backup power capacities.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
from _helpers import set_scenario_config
import cartopy.crs as ccrs
from pypsa.plot import add_legend_lines, add_legend_circles, add_legend_patches

BACKUP_TYPES = {
    'OCGT': "#e05b09",
    'OCGT methanol': "#a1ffe6",
    'H2 turbine': "#d5d6f5",
    'H2 Fuel Cell': "#a6a8ed",
    'urban central gas CHP': "#f7b7a3",
    'urban central gas CHP CC': "#e69487",
    'urban central solid biomass CHP': "#d5ca8d",
    'urban central solid biomass CHP CC': "#baa741",
    'battery discharger': "#ace37f",
    'home battery discharger': "#80c944",
}

NICE_NAMES = {
    'OCGT': "gas turbine",
    'OCGT methanol': "methanol turbine",
    'H2 turbine': "hydrogen turbine",
    'H2 Fuel Cell': "hydrogen CHP",
    'urban central gas CHP': "gas CHP",
    'urban central gas CHP CC': "gas CHP CC",
    'urban central solid biomass CHP': "biomass CHP",
    'urban central solid biomass CHP CC': "biomass CHP CC",
    'battery discharger': "utility-scale battery",
    'home battery discharger': "home battery",
}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_backup_power_map",
            opts="",
            clusters="115",
            ll="vopt",
            sector_opts="imp",
            planning_horizons="2050",
            configfiles="config/config.20240826-z1.yaml",
        )
    set_scenario_config(snakemake)

    plt.style.use(snakemake.input.rc)

    lw_factor = 6e3
    bs_factor = 3e4

    n = pypsa.Network(snakemake.input.network)

    bus_sizes = n.links.query("carrier in @BACKUP_TYPES.keys()").groupby([n.links.bus1.map(n.buses.location), n.links.carrier]).apply(lambda x: (x.p_nom_opt * x.efficiency).sum())

    n.mremove("Bus", n.buses.query("carrier != 'AC'").index)
    n.buses = n.buses.loc[n.buses.index.str.len() != 2]
    n.buses = n.buses.loc[~n.buses.index.str.contains("-")]

    prices = ((n.snapshot_weightings.generators @ n.buses_t.marginal_price) / n.snapshot_weightings.generators.sum())[n.buses.index]

    regions = gpd.read_file(snakemake.input.regions).set_index("name")

    proj = ccrs.EqualEarth()

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"projection": proj})
    regions.to_crs(proj.proj4_init).plot(
        ax=ax, column=prices, cmap='Purples',
        legend=True,
        vmin=30,
        vmax=90,
        legend_kwds={
            "label": r"time-averaged electricity price [€/MWh]",
            "shrink": 0.5,
            "pad": 0.02,
            "orientation": "vertical",
            "extend": "both",
        },
    )
    n.plot(
        ax=ax,
        margin=0.06,
        line_widths=n.lines.s_nom_opt / lw_factor,
        link_widths=n.links.p_nom_opt / lw_factor,
        line_colors='skyblue',
        link_colors='plum',
        geomap='50m',
        bus_sizes=bus_sizes / bs_factor,
        bus_colors=BACKUP_TYPES,
    )

    sizes = [10, 20]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / lw_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc=[0.625, 1.05],
        frameon=False,
        labelspacing=2,
        handletextpad=1,
        fontsize=11,
    )

    add_legend_lines(
        ax,
        sizes,
        labels,
        patch_kw=dict(color="plum", linestyle=":", gapcolor="skyblue"),
        legend_kw=legend_kw,
    )

    legend_kw = dict(
        loc=[0, 1],
        frameon=False,
        labelspacing=0.5,
        handletextpad=1,
        fontsize=11,
        ncol=2,
    )

    add_legend_patches(
        ax, list(BACKUP_TYPES.values()), list(NICE_NAMES.values()), legend_kw=legend_kw
    )

    sizes = [10, 30]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / bs_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc=[0.82, 1.05],
        frameon=False,
        labelspacing=2,
        handletextpad=1,
        fontsize=11,
    )

    add_legend_circles(
        ax, sizes, labels, patch_kw=dict(facecolor="lightgrey"), legend_kw=legend_kw
    )

    plt.savefig(snakemake.output[0], bbox_inches="tight")
