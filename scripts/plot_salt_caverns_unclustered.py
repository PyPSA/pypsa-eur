# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot unclustered salt caverns.
"""

import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_salt_caverns_unclustered",
            configfiles=["../../config/config.test.yaml"]
        )

    plt.style.use(snakemake.input.rc)

    caverns = gpd.read_file(snakemake.input.caverns)

    crs = ccrs.EqualEarth()

    caverns = caverns.to_crs(crs.proj4_init)

    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": crs})

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    caverns.plot(
        ax=ax,
        column="storage_type",
        cmap="tab10_r",
        legend=True,
        linewidth=0,
        legend_kwds=dict(
            title="Salt Caverns for\nHydrogen Storage", loc=(0.21, 0.82)
        ),
    )


    plt.xlim(-1e6, 2.6e6)
    plt.ylim(4.3e6, 7.8e6)

    ax.axis("off")
    ax.gridlines(linestyle=":")

    for fn in snakemake.output:
        plt.savefig(fn)
