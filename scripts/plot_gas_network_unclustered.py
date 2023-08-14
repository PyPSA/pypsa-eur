# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot unclustered gas transmission network.
"""

import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from build_gas_input_locations import build_gas_input_locations
from matplotlib.ticker import LogFormatter
from shapely import wkt

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("plot_gas_network_unclustered")

    plt.style.use(snakemake.input.rc)

    countries = snakemake.config["countries"]

    pts = build_gas_input_locations(
        snakemake.input.gem,
        snakemake.input.entry,
        snakemake.input.storage,
        countries,
    )

    sums = pts.groupby("type").capacity.sum()

    for t in ["lng", "production", "pipeline", "storage"]:
        pts.loc[pts["type"] == t, "capacity"] /= sums[t]

    pts["type"] = pts["type"].replace(
        dict(
            production="Fossil Extraction",
            lng="LNG Terminal",
            pipeline="Entrypoint",
            storage="Storage",
        )
    )

    df = pd.read_csv(snakemake.input.gas_network, index_col=0)
    for col in ["geometry"]:
        df[col] = df[col].apply(wkt.loads)

    df = gpd.GeoDataFrame(df, geometry="geometry", crs="EPSG:4326")

    crs = ccrs.EqualEarth()

    df = df.to_crs(crs.proj4_init)
    pts = pts.to_crs(crs.proj4_init)

    fig, ax = plt.subplots(figsize=(9.5, 9.5), subplot_kw={"projection": crs})

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    vmax = 100

    df.plot(
        ax=ax,
        column=df["p_nom"].div(1e3).fillna(0.0),
        linewidths=np.log(df.p_nom.div(df.p_nom.min())).fillna(0.0).div(3),
        cmap="Spectral_r",
        vmax=vmax,
        legend=True,
        legend_kwds=dict(label="Gas Pipeline Capacity [GW]", shrink=0.55, extend="max"),
        norm=mcolors.LogNorm(vmin=1, vmax=vmax),
    )

    pts.plot(
        ax=ax,
        column="type",
        markersize=pts["capacity"].fillna(0.0) * 750,
        alpha=0.8,
        legend=True,
        legend_kwds=dict(loc=[0.08, 0.86]),
        zorder=6,
    )

    ax.gridlines(linestyle=":")
    ax.axis("off")

    for fn in snakemake.output:
        plt.savefig(fn)
