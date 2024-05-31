# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates map of network with import options.
"""

import logging

logger = logging.getLogger(__name__)

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from _helpers import configure_logging
from plot_power_network import assign_location
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

PTX_THRESHOLD = 1  # M€


def import_costs(n):
    kwargs = dict(comps={"Generator"}, groupby=["bus", "carrier"])
    tsc_gen = n.statistics.capex(**kwargs) + n.statistics.opex(**kwargs)

    kwargs = dict(comps={"Link"}, groupby=["bus0", "carrier"])
    tsc_link = n.statistics.capex(**kwargs) + n.statistics.opex(**kwargs)
    tsc_link.rename_axis(index={"bus0": "bus"}, inplace=True)

    kwargs = dict(comps={"Store"}, groupby=["bus", "carrier"])
    tsc_sto = n.statistics.capex(**kwargs) + n.statistics.opex(**kwargs)

    # groupby needed to merge duplicate entries
    tsc_imp = (
        pd.concat([tsc_gen, tsc_link, tsc_sto])
        .droplevel(0)
        .groupby(["bus", "carrier"])
        .sum()
    ) / 1e6  # M€

    tsc_imp.index = pd.MultiIndex.from_arrays(
        [
            tsc_imp.index.get_level_values("bus").map(n.buses.location),
            tsc_imp.index.get_level_values("carrier"),
        ]
    )

    locs = tsc_imp.index.get_level_values("bus")
    carriers = tsc_imp.index.get_level_values("carrier")
    to_skip = ["", "EU", "EU import", "EU solid", "process"]
    import_components = carriers.str.contains(
        "(external|import shipping|import pipeline)"
    )
    tsc_imp = tsc_imp[~(locs.isna() | locs.isin(to_skip)) & import_components]
    tsc_imp.index = tsc_imp.index.remove_unused_levels()

    tsc_imp = tsc_imp.where(tsc_imp > PTX_THRESHOLD).dropna()

    return tsc_imp.sort_index()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_networks",
            simpl="",
            opts="",
            clusters="100",
            ll="vopt",
            sector_opts="Co2L0-2190SEG-T-H-B-I-S-A-nowasteheat",
            planning_horizons="2050",
            configfiles="../../config/config.20230926-zecm.yaml",
        )

    configure_logging(snakemake)

    plt.style.use(snakemake.input.rc)

    tech_colors = snakemake.config["plotting"]["tech_colors"]

    crs = ccrs.EqualEarth()

    bus_size_factor = 7e3
    line_width_factor = 1e4

    n = pypsa.Network(snakemake.input.network)
    assign_location(n)
    if n.links.carrier.str.contains("hvdc-to-elec").any():
        exporters = snakemake.config["sector"]["import"]["endogenous_hvdc_import"][
            "exporters"
        ]
        n.buses.loc[exporters, "location"] = exporters

    regions = (
        gpd.read_file(snakemake.input.regions, crs=4326)
        .set_index("name")
        .to_crs(crs.proj4_init)
    )

    ptx = n.links.filter(
        regex="(methanolisation|Haber|Fischer|Sabatier|Electrolysis)", axis=0
    ).copy()
    ptx["cost"] = ptx.eval("p_nom_opt * capital_cost / 1e6")
    ptx = ptx.groupby(ptx.bus0.map(n.buses.location)).cost.sum()
    ptx = ptx.loc[ptx > PTX_THRESHOLD]
    regions["ptx"] = ptx

    imp_costs = import_costs(n)

    # add EU import costs

    eu_sizes = {}

    carriers = ["EU import shipping-lnh3", "EU import shipping-steel", "EU import shipping-hbi"]
    for carrier in carriers:
        if carrier in n.generators.index:
            mc = n.generators.at[carrier, "marginal_cost"]
            energy = n.generators_t.p.loc[:, carrier] @ n.snapshot_weightings.generators
            eu_sizes[tuple(carrier.split(" ", 1))] = mc * energy / 1e6

    carriers = ["EU import shipping-ftfuel", "EU import shipping-meoh"]
    for carrier in carriers:
        if carrier in n.links.index:
            mc = n.links.at[carrier, "marginal_cost"]
            energy = n.links_t.p0.loc[:, carrier] @ n.snapshot_weightings.generators
            eu_sizes[tuple(carrier.split(" ", 1))] = mc * energy / 1e6

    # place EU bus in sea
    if eu_sizes:
        eu_sizes = pd.Series(eu_sizes)
        eu_sizes.index.names = ["bus", "carrier"]
        if "EU" in n.buses.index:
            n.remove("Bus", "EU")
        n.add("Bus", "EU", x=-9, y=47.5)

        imp_costs = pd.concat([imp_costs, eu_sizes])

    # patch network
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)
    if "KZ" in n.buses.index:
        n.buses.loc["KZ", "x"] = 52
        n.buses.loc["KZ", "y"] = 49
    if "CN-West" in n.buses.index:
        n.buses.loc["CN-West", "x"] = 79
        n.buses.loc["CN-West", "y"] = 38
    for ct in n.buses.index.intersection({"MA", "DZ", "TN", "LY", "EG", "SA"}):
        n.buses.loc[ct, "y"] += 2

    link_colors = pd.Series(
        n.links.index.map(
            lambda x: "olivedrab" if "import hvdc-to-elec" in x else "orange"
        ),
        index=n.links.index,
    )

    link_widths = (
        n.links.p_nom_opt.where(n.links.p_nom_opt > 1e3)
        .fillna(0.0)
        .div(line_width_factor)
    )
    line_widths = (
        n.lines.s_nom_opt.where(n.lines.s_nom_opt > 1e3)
        .fillna(0.0)
        .div(line_width_factor)
    )

    fig, ax = plt.subplots(subplot_kw={"projection": crs}, figsize=(8, 12))

    regions.plot(
        ax=ax,
        column="ptx",
        cmap="Blues",
        linewidths=0.5,
        edgecolor="lightgray",
        legend=True,
        vmin=0,
        vmax=6000,
        legend_kwds={
            "label": r"PtX investments [M€/a]",
            "shrink": 0.35,
            "pad": 0.015,
        },
    )

    n.plot(
        ax=ax,
        color_geomap={"ocean": "white", "land": "#efefef"},
        bus_sizes=imp_costs.div(bus_size_factor),
        bus_colors=tech_colors,
        line_colors="olive",
        line_widths=line_widths,
        link_widths=link_widths,
        link_colors=link_colors,
    )

    n.plot(
        ax=ax,
        bus_sizes=0.0,
        line_colors="tan",
        link_colors="tan",
        line_widths=n.lines.s_nom.div(line_width_factor),
        link_widths=n.links.p_nom.div(line_width_factor),
        boundaries=[-11, 48, 25.25, 71.5],
        margin=0,
    )

    patches = imp_costs.index.get_level_values(1).unique()
    labels = list(patches) + [
        "import HVDC (new)",
        "internal HVDC (new)",
        "internal HVAC (existing)",
        "internal HVAC (reinforced)",
    ]
    colors = [tech_colors[c] for c in patches] + ["olivedrab", "orange", "tan", "olive"]

    legend_kw = dict(
        bbox_to_anchor=(1.65, 0.6),
        frameon=False,
        title="technology",
        alignment="left",
    )

    add_legend_patches(
        ax,
        colors,
        labels,
        legend_kw=legend_kw,
    )

    legend_kw = dict(
        bbox_to_anchor=(1.52, 1.04),
        frameon=False,
        title="import expenditures",
        alignment="left",
    )

    add_legend_circles(
        ax,
        [1000 / bus_size_factor, 5000 / bus_size_factor, 10000 / bus_size_factor],
        ["1,000 M€/a", "5,000 M€/a", "10,000 M€/a"],
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    legend_kw = dict(
        bbox_to_anchor=(1.52, 0.83),
        frameon=False,
        title="transmission capacity",
        alignment="left",
    )

    add_legend_lines(
        ax,
        [
            10000 / line_width_factor,
            20000 / line_width_factor,
            30000 / line_width_factor,
        ],
        ["10 GW", "20 GW", "30 GW"],
        patch_kw=dict(color="lightgrey"),
        legend_kw=legend_kw,
    )

    for fn in snakemake.output:
        plt.savefig(fn, bbox_inches="tight")
