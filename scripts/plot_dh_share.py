# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot district heating shares.
"""

import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_dh_share",
            clusters=115,
            planning_horizons=2050,
            configfiles=["config/config.20240826-z1.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    dh_share = pd.read_csv(snakemake.input.dh_share, index_col=0)["district fraction of node"].mul(100)  # %

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    crs = ccrs.EqualEarth()

    regions = regions.to_crs(crs.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": crs})

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    regions.plot(
        ax=ax,
        column=dh_share.reindex(regions.index),
        cmap="Reds",
        legend=True,
        linewidth=0.5,
        vmax=60,
        edgecolor="grey",
        legend_kwds={
            "label": "district heating share of residential/services heat demand [%]",
            "shrink": 0.7,
            "extend": "max",
        },
    )

    plt.xlim(-1e6, 2.6e6)
    plt.ylim(4.3e6, 7.8e6)

    ax.axis("off")

    for fn in snakemake.output:
        plt.savefig(fn)