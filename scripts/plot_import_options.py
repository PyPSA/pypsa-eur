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
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from _helpers import configure_logging
from plot_power_network import assign_location
from pypsa.plot import add_legend_circles, add_legend_patches

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_options",
            simpl="",
            opts="",
            clusters="100",
            ll="v1.5",
            sector_opts="Co2L0-2190SEG-T-H-B-I-S-A-imp",
            planning_horizons="2050",
            configfiles="../../config/config.100n-seg.yaml",
        )

    configure_logging(snakemake)

    plt.style.use(snakemake.input.rc)

    tech_colors = snakemake.config["plotting"]["tech_colors"]
    tech_colors["lng"] = "tomato"
    tech_colors["pipeline"] = "orchid"

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
    regions["marginal_cost"] = h2_cost.groupby(
        h2_cost.bus.map(n.buses.location)
    ).marginal_cost.min()

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
            lambda x: "olivedrab" if "import hvdc-to-elec" in x else "tan"
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
        [pd.Series(0.3, index=mi), inputs.div(bus_size_factor)], axis=0
    )

    fig, ax = plt.subplots(subplot_kw={"projection": crs}, figsize=(8, 12))

    n.plot(
        ax=ax,
        color_geomap={"ocean": "white", "land": "#efefef"},
        bus_sizes=bus_sizes_plain,
        bus_colors=tech_colors,
        line_colors="tan",
        line_widths=0.75,
        link_widths=0.75,
        link_colors=link_colors,
        boundaries=[-11, 48, 26.5, 70],
        margin=0,
    )

    regions.plot(
        ax=ax,
        column="marginal_cost",
        cmap="Blues_r",
        linewidths=0,
        vmin=50,
        vmax=100,
        legend=True,
        legend_kwds={
            "label": r"H$_2$ import cost [€/MWh]",
            "shrink": 0.35,
            "pad": 0.015,
        },
    )

    names = {
        "external onwind": "onshore wind",
        "external solar-utility": "solar PV",
        "external battery": "battery storage",
        "external H2": "hydrogen storage",
    }
    labels = list(names.values()) + ["HVDC import link", "internal power line"]
    colors = [tech_colors[c] for c in names.keys()] + ["olivedrab", "tan"]

    legend_kw = dict(
        bbox_to_anchor=(1.51, 1.03), frameon=False, title="electricity imports"
    )

    add_legend_patches(
        ax,
        colors,
        labels,
        legend_kw=legend_kw,
    )

    legend_kw = dict(
        bbox_to_anchor=(1.46, 0.66), frameon=False, title="H$_2$ and CH$_4$ imports"
    )

    names = {
        "lng": "LNG terminal",
        "pipeline": "pipeline entry",
    }
    labels = list(names.values())
    colors = [tech_colors[c] for c in names.keys()]

    add_legend_patches(
        ax,
        colors,
        labels,
        legend_kw=legend_kw,
    )

    legend_kw = dict(
        bbox_to_anchor=(1.42, 0.47),
        frameon=False,
        title="import capacity",
    )

    add_legend_circles(
        ax,
        [100e3 / bus_size_factor],
        ["100 GW"],
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    cost_range = (
        pd.concat(
            [
                c.df.filter(like="import", axis=0)
                .groupby("carrier")
                .marginal_cost.describe()
                for c in n.iterate_components({"Link", "Generator"})
            ]
        )
        .drop("import hvdc-to-elec")
        .sort_values(by="min")[["min", "max"]]
        .astype(int)
        .T
    )

    translate = {
        "import shipping-ftfuel": "€/MWh FT",
        "import shipping-meoh": "€/MWh MeOh",
        "import pipeline-h2": r"€/MWh H2$_{(g)}$",
        "import shipping-lh2": r"€/MWh H2$_{(l)}$",
        "import shipping-lch4": r"€/MWh CH4$_{(l)}$",
        "import shipping-lnh3": r"€/MWh NH3$_{(l)}$",
        "import shipping-steel": r"€/t steel",
    }
    text = ""
    for carrier, values in cost_range.items():
        if abs(values["min"] - values["max"]) < 1:
            value = str(values["min"])
        else:
            value = str(values["min"]) + "-" + str(values["max"])
        text += value + " " + translate[carrier] + "\n"

    ax.text(1.2, 0.0, text, transform=ax.transAxes, linespacing=1.2)

    for fn in snakemake.output:
        plt.savefig(fn, bbox_inches="tight")
