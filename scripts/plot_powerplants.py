# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot power plants.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import matplotlib.colors as mcolors

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_powerplants",
            configfiles=["../../config/config.test.yaml"]
        )

    plt.style.use(snakemake.input.rc)

    df = pd.read_csv(snakemake.input.powerplants, index_col=0)

    df = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs="EPSG:4326")

    colors = {
        "Bioenergy": "#80c944",
        "Hard Coal": "black",
        "Hydro": "#235ebc",
        "Lignite": "#826837",
        "Natural Gas": "#a85522",
        "Nuclear": "#ff8c00",
        "Oil": "#c9c9c9",
        "Waste": "purple",
        "Other": "gold"
    }

    crs = ccrs.EqualEarth()

    df = df.cx[-12:30, 35:72]
    df = df.to_crs(crs.proj4_init)

    fig, ax = plt.subplots(figsize=(8,8), subplot_kw={"projection": crs})

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    df.plot(
        ax=ax,
        column="Fueltype",
        markersize=df["Capacity"] / 35,
        alpha=0.75,
        legend=True,
        cmap=mcolors.ListedColormap(
            pd.Series(df.Fueltype.unique()).sort_values().map(colors).values
        ),
        legend_kwds=dict(title="Technology (radius ~ capacity)", loc=[0.13, 0.82], ncols=2, title_fontproperties={'weight':'bold'}),
    )

    ax.axis("off")

    for fn in snakemake.output:
        plt.savefig(fn)
