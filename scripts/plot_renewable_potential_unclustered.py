# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot unclustered wind and solar renewable potentials.
"""

import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_renewable_potential_unclustered",
            configfiles=["../../config/config.test.yaml"],
        )

    plt.style.use(snakemake.input.rc)

    regions = pd.concat(
        [
            gpd.read_file(snakemake.input.regions_onshore),
            gpd.read_file(snakemake.input.regions_offshore),
        ]
    )
    regions = regions.dissolve("name")

    regions["Area"] = regions.to_crs(epsg=3035).area.div(1e6)

    wind = pd.Series()
    for profile in ["onwind", "offwind-ac", "offwind-dc"]:
        ds = xr.open_dataset(snakemake.input[f"profile_{profile}"])
        wind = pd.concat([wind, (ds.p_nom_max * ds.profile.sum("time")).to_pandas()])
    wind = wind.groupby(level=0).sum().reindex(regions.index, fill_value=0)
    wind_per_skm = wind / regions.Area / 1e3  # GWh

    proj = ccrs.EqualEarth()
    regions = regions.to_crs(proj.proj4_init)
    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": proj})
    regions.plot(
        ax=ax,
        column=wind_per_skm,
        cmap="Blues",
        linewidths=0,
        legend=True,
        legend_kwds={"label": r"Wind Energy Potential [GWh/a/km$^2$]", "shrink": 0.8},
    )
    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=4)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.gridlines(linestyle=":")
    ax.axis("off")

    for fn in snakemake.output["wind"]:
        plt.savefig(fn)

    onregions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
    onregions["Area"] = onregions.to_crs(epsg=3035).area.div(1e6)

    ds = xr.open_dataset(snakemake.input.profile_solar)
    solar = (ds.p_nom_max * ds.profile.sum("time")).to_pandas()

    solar = solar.groupby(level=0).sum().reindex(onregions.index, fill_value=0)
    solar_per_skm = solar / onregions.Area / 1e3  # GWh

    proj = ccrs.EqualEarth()
    onregions = onregions.to_crs(proj.proj4_init)
    fig, ax = plt.subplots(figsize=(7, 7), subplot_kw={"projection": proj})
    onregions.plot(
        ax=ax,
        column=solar_per_skm,
        cmap="Reds",
        linewidths=0,
        legend=True,
        legend_kwds={"label": r"Solar Energy Potential [GWh/a/km$^2$]", "shrink": 0.8},
    )
    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.gridlines(linestyle=":")
    ax.axis("off")

    for fn in snakemake.output["solar"]:
        plt.savefig(fn)
