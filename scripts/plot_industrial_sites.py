# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Plot industrial sites.
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import country_converter as coco
cc = coco.CountryConverter()


def prepare_hotmaps_database():
    """
    Load hotmaps database of industrial sites.
    """

    df = pd.read_csv(
        snakemake.input.hotmaps, sep=";", index_col=0
    )

    df[["srid", "coordinates"]] = df.geom.str.split(";", expand=True)

    # remove those sites without valid locations
    df.drop(df.index[df.coordinates.isna()], inplace=True)

    df["coordinates"] = gpd.GeoSeries.from_wkt(df["coordinates"])
    
    gdf = gpd.GeoDataFrame(df, geometry="coordinates", crs="EPSG:4326")


    return gdf


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_industrial_sites",
            configfiles=["../../config/config.test.yaml"]
        )

    plt.style.use(snakemake.input.rc)

    crs = ccrs.EqualEarth()

    countries = gpd.read_file(snakemake.input.countries).set_index('name')

    hotmaps = prepare_hotmaps_database()
    hotmaps = hotmaps.cx[-12:30, 35:72]
    hotmaps = hotmaps.to_crs(crs.proj4_init)

    not_represented = ["AL", "BA", "RS", "MK", "ME"]
    missing_countries = countries.loc[countries.index.intersection(not_represented)]

    fig, ax = plt.subplots(figsize=(8,8), subplot_kw={"projection": crs})

    ax.add_feature(cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS.with_scale("50m"), linewidth=0.5, zorder=2)

    missing_countries.to_crs(crs.proj4_init).plot(
        ax=ax,
        color='lightgrey',
    )

    emissions = hotmaps["Emissions_ETS_2014"].fillna(hotmaps["Emissions_EPRTR_2014"]).fillna(hotmaps["Production"]).fillna(2e4)

    hotmaps.plot(
        ax=ax,
        column="Subsector",
        markersize=emissions / 4e4,
        alpha=0.5,
        legend=True,
        legend_kwds=dict(title="Industry Sector (radius ~ emissions)", frameon=False, ncols=2, loc=[0,.85], title_fontproperties={'weight':'bold'}),
    )
    ax.axis("off")

    for fn in snakemake.output:
        plt.savefig(fn)
