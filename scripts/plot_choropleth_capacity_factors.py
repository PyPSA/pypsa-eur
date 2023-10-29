# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot average capacity factor map.
"""

import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pypsa
from _helpers import ensure_output_dir_exists


def plot_choropleth(
    df,
    geodf,
    carrier,
    cmap="Blues",
    vmin=0,
    vmax=100,
    label="capacity factors [%]",
    dir=".",
    legend_kwds=None,
    title="",
    **kwargs,
):
    if legend_kwds is None:
        legend_kwds = {
            "label": label,
            "shrink": 0.7,
            "extend": "max",
        }

    proj = ccrs.EqualEarth()
    geodf = geodf.to_crs(proj.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": proj})

    geodf.plot(
        ax=ax,
        column=df[carrier].reindex(geodf.index),
        cmap=cmap,
        linewidths=0.5,
        legend=True,
        vmax=vmax,
        vmin=vmin,
        edgecolor="lightgray",
        legend_kwds=legend_kwds,
        **kwargs,
    )

    ax.set_title(title)

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    ax.axis("off")

    carrier_fn = carrier.replace("-", "_").replace(" ", "_")
    fn = f"map-{carrier_fn}"
    if not dir.endswith("-"):
        dir += "/"
    plt.savefig(dir + fn + ".pdf", bbox_inches="tight")
    plt.savefig(dir + fn + ".png", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_choropleth_capacity_factors",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    ensure_output_dir_exists(snakemake)
    dir = snakemake.output[0]

    n = pypsa.Network(snakemake.input.network)

    df = (
        n.generators_t.p_max_pu.mean()
        .groupby([n.generators.carrier, n.generators.bus])
        .first()
        .unstack(0)
        .mul(100)
    )

    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    regions_offshore = gpd.read_file(snakemake.input.regions_offshore).set_index("name")

    plot_choropleth(df, regions_onshore, "onwind", vmax=50, dir=dir)
    plot_choropleth(df, regions_onshore, "solar", "Oranges", vmax=15, dir=dir)
    plot_choropleth(df, regions_onshore, "ror", "GnBu", dir=dir)

    plot_choropleth(df, regions_offshore, "offwind-dc", vmax=50, dir=dir)
    plot_choropleth(df, regions_offshore, "offwind-ac", vmax=50, dir=dir)
