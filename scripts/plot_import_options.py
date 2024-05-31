# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Creates map of import options.
"""

import logging

logger = logging.getLogger(__name__)

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import country_converter as coco
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
import seaborn as sns
from _helpers import configure_logging
from plot_power_network import assign_location
from pypsa.plot import add_legend_circles, add_legend_patches

cc = coco.CountryConverter()

# for EU: https://ec.europa.eu/eurostat/databrowser/view/prc_hicp_aind__custom_9900786/default/table?lang=en
EUR_2015_TO_2020 = 1.002 * 1.017 * 1.019 * 1.015 * 1.007

NICE_NAMES = {
    "pipeline-h2": r"H$_2$ (pipeline)",
    "shipping-lh2": "H$_2$ (ship)",
    "shipping-lnh3": "ammonia",
    "shipping-lch4": "methane",
    "shipping-meoh": "methanol",
    "shipping-ftfuel": "Fischer-Tropsch",
    "shipping-steel": "steel",
    "shipping-hbi": "hbi",
}

PALETTE = {
    "Argentina": "#74acdf",
    "Algeria": "#d21034",
    "Namibia": "#003580",
    "Saudi Arabia": "#006c35",
    "Chile": "darkorange",
    "Other": "#aaa",
}


def create_stripplot(ic, ax):

    order = list(NICE_NAMES.values())[:-1]
    minimums = (
        ic.groupby("esc").value.min().round(1)[order].reset_index(drop=True).to_dict()
    )
    maximums = (
        ic.groupby("esc").value.max().round(1)[order].reset_index(drop=True).to_dict()
    )

    sns.stripplot(
        data=ic,
        x="esc",
        y="value",
        alpha=0.6,
        hue="exporter",
        jitter=0.28,
        palette=PALETTE,
        ax=ax,
        order=order,
        size=4,
    )
    sns.violinplot(
        data=ic,
        x="esc",
        y="value",
        linewidth=0,
        saturation=0.3,
        cut=0,
        color="#ddd",
        fill=True,
        ax=ax,
        order=order,
        zorder=-1,
    )

    ax.set_ylim(0, 200)
    ax.set_xlabel("")
    ax.set_ylabel("import cost [€/MWh]")
    ax.grid(False)
    ax.set_yticks(np.arange(0, 210, 20))
    ax.set_yticks(np.arange(10, 210, 20), minor=True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=18, ha="right")
    handles, labels = ax.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()
    new_handles = []
    for handle in handles:
        handle.set_markersize(handle.get_markersize() * 2)
        new_handles.append(handle)
    handles = new_handles
    for x, y in minimums.items():
        ax.text(x, y - 10, str(y), ha="center", va="bottom", fontsize=9)
    for x, y in maximums.items():
        ax.text(x, y + 5, str(y), ha="center", va="bottom", fontsize=9)
    ax.legend(title="", ncol=1, loc=(0.55, 0.05), labelspacing=0.3, frameon=False)
    for spine in ax.spines.values():
        spine.set_visible(False)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_options",
            simpl="",
            opts="",
            clusters="110",
            ll="vopt",
            sector_opts="Co2L0-2190SEG-T-H-B-I-S-A-imp+AC",
            planning_horizons="2050",
            configfiles="../../config/config.20231025-zecm.yaml",
        )

    configure_logging(snakemake)

    plt.style.use(["bmh", snakemake.input.rc])

    # dummy output if no imports considered
    if "-imp" not in snakemake.wildcards.sector_opts:
        import sys

        fig, ax = plt.subplots()
        for fn in snakemake.output:
            plt.savefig(fn, bbox_inches="tight")
        sys.exit(0)

    tech_colors = snakemake.config["plotting"]["tech_colors"]
    tech_colors["lng"] = "#e37959"
    tech_colors["pipeline"] = "#86cfbc"

    crs = ccrs.EqualEarth()

    bus_size_factor = 7.5e4

    n = pypsa.Network(snakemake.input.network)
    assign_location(n)

    regions = (
        gpd.read_file(snakemake.input.regions).set_index("name").to_crs(crs.proj4_init)
    )

    inputs = pd.read_csv(snakemake.input.entrypoints, index_col=0)[
        ["lng", "pipeline"]
    ].copy()
    countries = ["DE", "GB", "BE", "FR", "EE", "LV", "LT", "FI"]
    pattern = "|".join(countries)
    inputs.loc[inputs.index.str.contains(pattern), "pipeline"] = 0.0
    inputs = inputs.stack()

    # TODO size external nodes according to wind and solar potential

    h2_cost = n.generators.filter(regex="import (pipeline-h2|shipping-lh2)", axis=0)
    regions["marginal_cost"] = (
        h2_cost.groupby(h2_cost.bus.map(n.buses.location)).marginal_cost.min()
        * EUR_2015_TO_2020
    )

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
            lambda x: "seagreen" if "import hvdc-to-elec" in x else "#b18ee6"
        ),
        index=n.links.index,
    )

    exporters = snakemake.config["sector"]["import"]["endogenous_hvdc_import"][
        "exporters"
    ]
    techs = [
        "external solar-utility",
        "external onwind",
        "external battery",
        "external H2",
    ]
    mi = pd.MultiIndex.from_product([exporters, techs])
    bus_sizes_plain = pd.concat(
        [pd.Series(0.4, index=mi), inputs.div(bus_size_factor)], axis=0
    )

    df = pd.read_csv(snakemake.input.imports, sep=";", keep_default_na=False)

    df["exporter"] = df.exporter.replace("", "NA")
    ic = df.query("subcategory == 'Cost per MWh delivered' and esc != 'hvdc-to-elec'")
    ic["exporter"] = ic.exporter.str.split("-").str[0]

    highlighted_countries = ["DZ", "AR", "SA", "CL"]
    ic["exporter"] = ic.exporter.apply(
        lambda x: (
            cc.convert(names=x, to="name_short")
            if x in highlighted_countries
            else "Other"
        )
    )

    ic["esc"] = ic.esc.map(NICE_NAMES)
    ic["value"] *= EUR_2015_TO_2020

    fig, ax = plt.subplots(subplot_kw={"projection": crs}, figsize=(10.5, 14))

    n.plot(
        ax=ax,
        color_geomap={"ocean": "white", "land": "#efefef"},
        bus_sizes=bus_sizes_plain,
        bus_colors=tech_colors,
        line_colors="#b18ee6",
        line_widths=1,
        link_widths=1,
        link_colors=link_colors,
        boundaries=[-11, 48, 25.25, 71.5],
        margin=0,
    )

    regions.plot(
        ax=ax,
        column="marginal_cost",
        cmap="Blues_r",
        edgecolor="#ddd",
        linewidths=0.5,
        vmin=50,
        vmax=100,
        legend=True,
        legend_kwds={
            "label": r"H$_2$ import cost [€/MWh]",
            "shrink": 0.53,
            "pad": 0.015,
            "aspect": 35,
        },
    )

    names = {
        "external onwind": "onshore wind",
        "external solar-utility": "solar PV",
        "external battery": "battery storage",
        "external H2": "hydrogen storage",
    }
    labels = list(names.values()) + [
        "HVDC import link",
        "internal power line",
        "LNG terminal",
        "pipeline entry",
    ]
    colors = [tech_colors[c] for c in names.keys()] + [
        "seagreen",
        "#b18ee6",
        "#e37959",
        "#86cfbc",
    ]

    legend_kw = dict(
        loc=(0.595, 0.87),
        frameon=True,
        title="",
        ncol=2,
        framealpha=1,
        borderpad=0.5,
        facecolor="white",
    )

    add_legend_patches(
        ax,
        colors,
        labels,
        legend_kw=legend_kw,
    )

    legend_kw = dict(
        loc=(0.623, 0.775),
        frameon=True,
        title="existing gas import capacity",
        ncol=3,
        labelspacing=1.1,
        framealpha=1,
        borderpad=0.5,
        facecolor="white",
    )

    add_legend_circles(
        ax,
        [10e3 / bus_size_factor, 50e3 / bus_size_factor, 100e3 / bus_size_factor],
        ["10 GW", "50 GW", "100 GW"],
        patch_kw=dict(facecolor="#86cfbc"),
        legend_kw=legend_kw,
    )

    ax.add_feature(
        cfeature.BORDERS.with_scale("50m"),
        linewidth=0.75,
        color="k",
    )

    ax.add_feature(
        cfeature.COASTLINE.with_scale("50m"),
        linewidth=0.75,
        color="k",
    )

    plt.tight_layout()

    for fn in snakemake.output.map:
        plt.savefig(fn, bbox_inches="tight")

    fig, ax_lr = plt.subplots(figsize=(3.9, 7.2))
    create_stripplot(ic, ax_lr)
    plt.tight_layout()
    for fn in snakemake.output.distribution:
        plt.savefig(fn, bbox_inches="tight")
