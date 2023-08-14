# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot unclustered salt caverns.
"""

import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd


def plot_biomass_potentials(bio, regions, kind):
    crs = ccrs.EqualEarth()
    regions = regions.to_crs(crs.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": crs})

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    nkind = "disregarded biomass" if kind == "not included" else kind
    label = f"{nkind} potentials [TWh/a]"

    regions.plot(
        ax=ax,
        column=bio[kind],
        cmap="Greens",
        legend=True,
        linewidth=0.5,
        edgecolor="grey",
        legend_kwds={
            "label": label,
            "shrink": 0.7,
            "extend": "max",
        },
    )

    total = bio[kind].sum()

    ax.text(-0.8e6, 7.4e6, f"{total:.0f} TWh/a", color="#343434")

    plt.xlim(-1e6, 2.6e6)
    plt.ylim(4.3e6, 7.8e6)

    ax.axis("off")

    for fn in snakemake.output[kind.replace(" ", "_")]:
        plt.savefig(fn)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_biomass_potentials",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    bio = pd.read_csv(snakemake.input.biomass, index_col=0).div(1e6)  # TWh/a

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    plot_biomass_potentials(bio, regions, kind="not included")
    plot_biomass_potentials(bio, regions, kind="biogas")
    plot_biomass_potentials(bio, regions, kind="solid biomass")
