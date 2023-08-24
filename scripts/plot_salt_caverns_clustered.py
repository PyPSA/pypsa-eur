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
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd


def plot_salt_caverns_by_node(
    cavern_nodes,
    cavern_regions,
    storage_type="onshore",
    cmap="GnBu",
    vmin=1,
    vmax=3000,
    label=r"H$_2$ Storage Potential [TWh]",
):
    crs = ccrs.EqualEarth()

    cavern_regions = cavern_regions.to_crs(crs.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": crs})

    cavern_regions.plot(
        ax=ax,
        column=cavern_nodes[storage_type].reindex(cavern_regions.index),
        cmap=cmap,
        legend=True,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        edgecolor="darkgray",
        legend_kwds={
            "label": label,
            "shrink": 0.7,
            "extend": "max",
        },
        norm=mcolors.LogNorm(vmin=1, vmax=vmax),
    )

    plt.title(f"{storage_type.capitalize()} Salt Cavern H$_2$ Storage Potentials")

    plt.xlim(-1.2e6, 2.6e6)
    plt.ylim(4.3e6, 7.8e6)

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    ax.axis("off")

    for fn in snakemake.output[storage_type]:
        plt.savefig(fn)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_salt_caverns_clustered",
            clusters=128,
            configfiles=["../../config/config.test.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    cavern_nodes = pd.read_csv(snakemake.input.caverns, index_col=0)
    cavern_nodes = cavern_nodes.where(cavern_nodes > 0.5)

    cavern_regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    cavern_offregions = gpd.read_file(snakemake.input.regions_offshore).set_index(
        "name"
    )

    plot_salt_caverns_by_node(cavern_nodes, cavern_regions, storage_type="onshore")
    plot_salt_caverns_by_node(cavern_nodes, cavern_regions, storage_type="nearshore")
    plot_salt_caverns_by_node(cavern_nodes, cavern_offregions, storage_type="offshore")
