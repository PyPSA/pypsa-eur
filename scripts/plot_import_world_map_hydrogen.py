# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023-2024 Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Creates world map of global green hydrogen export costs.
"""

import logging

import cartopy.crs as ccrs
import country_converter as coco
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from _helpers import configure_logging
from shapely.geometry import box


from plot_import_world_map import add_legend_patches, reduce_resolution

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()

AREA_CRS = "ESRI:54009"


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_world_map_hydrogen",
            opts="",
            clusters=115,
            ll="vopt",
            sector_opts="imp",
            planning_horizons="2050",
            configfiles="config/config.20240826-z1.yaml",
        )

    configure_logging(snakemake)

    import_fn = snakemake.input.imports
    countries_fn = snakemake.input.countries

    config = snakemake.config

    plt.style.use(["bmh", snakemake.input.rc])

    # download world country shapes, version clipped to Europe, and GADM in AR

    world = gpd.read_file(countries_fn).set_index("ADM0_A3_DE")
    world.drop("ATA", inplace=True)
    world.rename(index={"KOS": "XKX"}, inplace=True)

    eu_countries = cc.convert(config["countries"], src="iso2", to="iso3")
    europe = world.loc[eu_countries]

    bounding_box = box(-12, 33, 42, 72)
    europe = gpd.clip(europe, bounding_box)

    # load import costs

    df = (
        pd.read_parquet(import_fn)
        .reset_index()
        .query("scenario == 'default' and year == 2040")
    )

    h2_vectors = ["shipping-lh2", "pipeline-h2"]
    import_costs = (
        df.query("subcategory == 'Cost per MWh delivered' and esc in @h2_vectors")
        .groupby(["exporter", "esc"])
        .value.min()
        .unstack()
    )

    import_vectors = import_costs.idxmin(axis=1)
    import_costs = import_costs.min(axis=1)

    import_regions = pd.concat(
        {
            idx: gpd.read_file(f"../data/imports/regions/{idx}.gpkg")
            .set_index("index")
            .loc["onshore"]
            for idx in df.exporter.unique()
        }
    )
    import_regions.index = import_regions.index.get_level_values(0)
    import_regions = gpd.GeoDataFrame(geometry=import_regions, crs=4326)
    import_regions = reduce_resolution(import_regions, tolerance=0.1, min_area=0.1)

    # create plot

    crs = ccrs.EqualEarth()

    fig, ax = plt.subplots(figsize=(14, 14), subplot_kw={"projection": crs})

    # main axis: choropleth layer

    world.to_crs(crs).plot(
        linewidth=1,
        edgecolor="black",
        ax=ax,
        color="#eee",
    )

    europe.to_crs(crs).plot(
        linewidth=1,
        edgecolor="black",
        ax=ax,
        color="#cbc7f0",
    )

    pipe_i = import_vectors.loc[import_vectors == "pipeline-h2"].index
    ship_i = import_vectors.loc[import_vectors != "pipeline-h2"].index

    import_regions.loc[pipe_i].to_crs(crs).plot(
        column=import_costs.reindex(import_regions.loc[pipe_i].index),
        hatch="..",
        linewidth=1,
        edgecolor="black",
        ax=ax,
        cmap="viridis_r",
        vmin=75,
        vmax=145,
    )

    import_regions.loc[ship_i].to_crs(crs).plot(
        column=import_costs.reindex(import_regions.loc[ship_i].index),
        linewidth=1,
        edgecolor="black",
        ax=ax,
        cmap="viridis_r",
        legend=True,
        vmin=75,
        vmax=145,
        legend_kwds=dict(
            label="Cost of hydrogen delivered to Europe [€/MWh]",
            orientation="horizontal",
            extend="both",
            shrink=0.6,
            aspect=30,
            pad=0.01,
        ),
    )

    add_legend_patches(
        ax,
        ["#eee", "#28b594", "#cbc7f0"],
        ["region not considered for export", "hydrogen export by pipeline", "region in European model scope"],
        ["", "..", ""],
        legend_kw=dict(
            bbox_to_anchor=(1, 0),
            frameon=False,
        ),
        patch_kw=dict(alpha=0.99),
    )

    import_regions.representative_point().to_crs(crs).annotate()

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_facecolor("none")
    fig.set_facecolor("none")

    ax.text(
        0.93,
        0.01,
        "Projection:\nEqual Earth",
        transform=ax.transAxes,
        fontsize=9,
        color="grey",
    )

    plt.tight_layout()

    for fn in snakemake.output:
        plt.savefig(fn, bbox_inches="tight")
