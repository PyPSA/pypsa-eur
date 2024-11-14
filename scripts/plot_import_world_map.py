# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- Fabian Neumann
#
# SPDX-License-Identifier: MIT
"""
Creates world map of global green fuel and material markets.
"""

import logging
import os

import cartopy.crs as ccrs
import country_converter as coco
import geopandas as gpd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import configure_logging
from atlite.gis import ExclusionContainer, shape_availability
from matplotlib.patches import Patch
from rasterio.features import geometry_mask
from rasterio.plot import show
from shapely.geometry import box


def add_legend_patches(ax, colors, labels, hatches, patch_kw=None, legend_kw=None):
    colors = np.atleast_1d(colors)
    labels = np.atleast_1d(labels)
    if patch_kw is None:
        patch_kw = {}
    if legend_kw is None:
        legend_kw = {}

    if len(colors) != len(labels):
        msg = "Colors and labels must have the same length."
        raise ValueError(msg)

    handles = [Patch(facecolor=c, hatch=h, **patch_kw) for c, h in zip(colors, hatches)]

    legend = ax.legend(handles, labels, **legend_kw)

    ax.get_figure().add_artist(legend)

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()

AREA_CRS = "ESRI:54009"

ARROW_COLOR = "#A0A0A0"

NICE_NAMES = {
    "pipeline-h2": r"H$_2$ (pipeline)",
    "shipping-lh2": "H$_2$ (ship)",
    "shipping-ftfuel": "Fischer-Tropsch",
    "shipping-meoh": "methanol",
    "shipping-lch4": "methane",
    "shipping-lnh3": "ammonia",
    "shipping-steel": "steel",
    "shipping-hbi": "HBI",
}


def reduce_resolution(gdf, tolerance=0.01, min_area=0.1):
    def simplify_and_filter(geom):
        simplified = geom.simplify(tolerance)
        if simplified.geom_type == 'MultiPolygon':
            return gpd.GeoSeries([poly for poly in simplified.geoms if poly.area > min_area]).union_all()
        return simplified if simplified.area > min_area else None

    # Create a copy of the input GeoDataFrame
    gdf_reduced = gdf.copy()
    
    # Apply the simplification and filtering
    gdf_reduced['geometry'] = gdf_reduced['geometry'].apply(simplify_and_filter)
    
    # Remove any rows with null geometries
    gdf_reduced = gdf_reduced.dropna(subset=['geometry'])
    
    return gdf_reduced


def rename(s):
    if "solar" in s:
        return "solar"
    if "wind" in s:
        return "wind"
    if "inverter" in s or "battery" in s:
        return "battery"
    if "storage" in s or "Buffer" in s or "storing" in s:
        return "fuel storage"
    if "transport" in s or "shipping fuel" in s or "dry bulk" in s or "pipeline" in s or "ship" in s or "HVDC" in s:
        return "transport"
    if "evaporation" in s or "liquefaction" in s or "compress" in s:
        return "evaporation/liquefaction"
    if "direct air capture" in s or "heat pump" in s:
        return "direct air capture"
    if s in [
        "Haber-Bosch (exp)",
        "Fischer-Tropsch (exp)",
        "methanolisation (exp)",
        "methanation (exp)",
    ]:
        return "hydrogen conversion"
    if "iron ore" in s:
        return "iron ore"
    if "direct iron reduction" in s:
        return "direct iron reduction"
    if "electric arc furnace" in s:
        return "electric arc furnace"
    return s.replace(" (exp)", "")


def get_cost_composition(df, country, escs, production):
    query_str = "category == 'cost' and exporter == @country and esc in @escs"
    composition = (
        df.query(query_str).groupby(["esc", "subcategory", "importer"]).value.min()
    )

    minimal = {}
    for name, group in composition.groupby("esc"):
        c = group.unstack("importer").droplevel("esc")
        minimal[name] = c[c.sum().idxmin()]
    composition = pd.concat(minimal, axis=1)

    composition = composition.groupby(rename).sum().div(production)

    composition = composition.where(composition > 0.01).dropna(how="all")

    sort_by = composition.sum().sort_values(ascending=True).index
    selection = pd.Index(COLORS.keys()).intersection(composition.index)
    composition = composition.loc[selection, sort_by].rename(columns=NICE_NAMES)
    return composition


def add_land_eligibility_example(ax, shape, glc_fn, wdpa_fn):

    excluder = ExclusionContainer(crs=AREA_CRS, res=200)
    excluder.add_raster(glc_fn, codes=[20, 30, 40, 60, 100], invert=True)
    excluder.add_raster(glc_fn, codes=[50], buffer=1500)
    if os.path.exists(wdpa_fn):
        wdpa = gpd.read_file(
            wdpa_fn,
            bbox=shape.to_crs(4326).geometry,
            layer="WDPA_Mar2024_Public_shp-polygons",
        ).to_crs(AREA_CRS)
        excluder.add_geometry(wdpa.geometry, buffer=1000)

    band, transform = shape_availability(shape, excluder)
    mask = ~geometry_mask(
        [shape.geometry.values[0]],
        transform=transform,
        invert=False,
        out_shape=band.shape,
    )
    masked_band = np.where(mask, ~band, np.nan)

    shape.plot(ax=ax, color="none", edgecolor="k", linewidth=1)
    show(masked_band, transform=transform, cmap="Purples", ax=ax)
    ax.set_axis_off()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_world_map",
            opts="",
            clusters=115,
            ll="vopt",
            sector_opts="imp",
            planning_horizons="2050",
            configfiles="config/config.20240826-z1.yaml",
        )

    configure_logging(snakemake)

    import_fn = snakemake.input.imports
    profile_fn = snakemake.input.profiles
    gadm_fn = snakemake.input.gadm_arg
    glc_fn = snakemake.input.copernicus_glc
    wdpa_fn = snakemake.input.wdpa
    countries_fn = snakemake.input.countries

    config = snakemake.config

    plt.style.use(["bmh", snakemake.input.rc])

    tech_colors = config["plotting"]["tech_colors"]

    COLORS = {
        "wind": tech_colors["onwind"],
        "solar": tech_colors["solar"],
        "battery": tech_colors["battery"],
        "electrolysis": tech_colors["H2 Electrolysis"],
        "fuel storage": '#ffd4dc',
        "hydrogen conversion": tech_colors["Fischer-Tropsch"],
        "direct air capture": tech_colors["DAC"],
        "iron ore": "#4e4f55",
        "direct iron reduction": tech_colors["steel"],
        "electric arc furnace": "#8795a8",
        "evaporation/liquefaction": "#8487e8",
        "transport": "#e0ae75",
    }

    # load capacity factor time series

    ds = xr.open_dataset(profile_fn)
    profile = ds.sel(exporter="MA", importer="EUSW").p_max_pu.to_pandas().T
    profile.rename(columns={"onwind": "wind", "solar-utility": "solar"}, inplace=True)

    # download world country shapes, version clipped to Europe, and GADM in AR

    world = gpd.read_file(countries_fn).set_index("ADM0_A3_DE")
    world.drop("ATA", inplace=True)
    world.rename(index={"KOS": "XKX"}, inplace=True)

    eu_countries = cc.convert(config["countries"], src="iso2", to="iso3")
    europe = world.loc[eu_countries]

    bounding_box = box(-12, 33, 42, 72)
    europe = gpd.clip(europe, bounding_box)

    gadm = gpd.read_file(gadm_fn, layer="ADM_ADM_1").set_index("NAME_1")
    shape = gadm.to_crs(AREA_CRS).loc[["Buenos Aires"]].geometry

    # load import costs

    df = (
        pd.read_parquet(import_fn)
        .reset_index()
        .query("scenario == 'default' and year == 2040")
    )

    import_costs = (
        df.query("subcategory == 'Cost per MWh delivered' and esc == 'shipping-meoh'")
        .groupby("exporter")
        .value.min()
    )

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

    composition_arg = get_cost_composition(
        df,
        "AR-South",
        [
            "shipping-lh2",
            "shipping-ftfuel",
            "shipping-meoh",
            "shipping-lch4",
            "shipping-lnh3",
        ],
        500e6,
    )

    composition_sau = get_cost_composition(
        df, "SA", ["pipeline-h2", "shipping-lh2"], 500e6
    )

    composition_aus = get_cost_composition(
        df, "AU-West", ["shipping-hbi", "shipping-steel"], 100e6
    )

    # add to legend

    composition_arg.loc["iron ore"] = pd.NA
    composition_arg.loc["direct iron reduction"] = pd.NA
    composition_arg.loc["electric arc furnace"] = pd.NA

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

    import_regions.to_crs(crs).plot(
        column=import_costs.reindex(import_regions.index),
        linewidth=1,
        edgecolor="black",
        ax=ax,
        cmap="viridis_r",
        legend=True,
        vmin=120,
        vmax=170,
        legend_kwds=dict(
            label="Cost for methanol fuel delivered to Europe [€/MWh]",
            orientation="horizontal",
            extend="max",
            shrink=0.6,
            aspect=30,
            pad=0.01,
        ),
        missing_kwds=dict(color="#ccc", hatch=".."),
    )

    add_legend_patches(
        ax,
        ["#eee", "#ccc", "#cbc7f0"],
        ["region not considered for export", "region with excluded ship exports", "region in European model scope"],
        ["", "..", ""],
        legend_kw=dict(
            bbox_to_anchor=(1, 0),
            frameon=False,
        ),
        patch_kw=dict(alpha=0.99)
    )

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

    # inset: wind and solar profiles

    ax_prof = ax.inset_axes([0.36, 0.68, 0.09, 0.11])

    week_profile = profile.loc["2013-03-01":"2013-03-07", ["solar", "wind"]]

    week_profile.plot(
        ax=ax_prof,
        linewidth=1,
        color=["gold", "royalblue"],
        ylim=(0, 1),
        clip_on=False,
    )

    ax_prof.legend(title="", loc=(0, 1), fontsize=8, ncol=2, columnspacing=0.8)
    ax_prof.set_xlabel("Day of March 2013", fontsize=8)
    ax_prof.set_ylabel("profile [p.u.]", fontsize=8)
    ax_prof.tick_params(axis="both", labelsize=8)

    ax_prof.xaxis.set_major_locator(mdates.DayLocator())
    ax_prof.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    xticks = week_profile.resample("D").mean().index
    ax_prof.set_xticks(xticks)
    ax_prof.set_xticklabels(xticks.day)

    for spine in ax_prof.spines.values():
        spine.set_visible(False)

    for label in ax_prof.get_xticklabels():
        label.set_fontsize(8)

    ax.annotate(
        "",
        xy=(0.45, 0.75),
        xytext=(0.485, 0.72),
        xycoords="axes fraction",
        arrowprops=dict(
            edgecolor=ARROW_COLOR,
            facecolor=ARROW_COLOR,
            linewidth=1.5,
            arrowstyle="-|>",
            connectionstyle="arc3,rad=0.2",
        ),
    )

    # inset: Argentina e-fuel import costs

    ax_arg = ax.inset_axes([0.07, 0.08, 0.09, 0.33])

    composition_arg.T.plot.bar(ax=ax_arg, stacked=True, color=COLORS)

    handles, labels = ax_arg.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()

    ax_arg.legend(handles, labels, title="", ncol=1, fontsize=9, loc=(1, 0))

    ax_arg.set_title("Import costs from\nSouthern Argentina\nto Europe by ship", fontsize=9)

    ax_arg.set_xlabel("")
    ax_arg.set_ylim(0, 140)
    ax_arg.set_yticks(range(0, 141, 20))
    ax_arg.set_yticks(range(10, 141, 20), minor=True)
    ax_arg.set_ylabel("€/MWh", fontsize=10)
    ax_arg.grid(axis="x")
    for spine in ax_arg.spines.values():
        spine.set_visible(False)

    ax.annotate(
        "",
        xy=(0.25, 0.1),
        xytext=(0.33, 0.15),
        xycoords="axes fraction",
        arrowprops=dict(
            edgecolor=ARROW_COLOR,
            facecolor=ARROW_COLOR,
            linewidth=1.5,
            arrowstyle="-|>",
            connectionstyle="arc3,rad=-0.2",
        ),
    )

    # inset: Saudi Arabia hydrogen pipeline versus ship imports

    ax_sau = ax.inset_axes([0.6625, 0.33, 0.036, 0.2])

    composition_sau.T.plot.bar(ax=ax_sau, stacked=True, color=COLORS, legend=False)

    ax_sau.set_title(r"LH$_2$" + " ship\nvs. pipeline", fontsize=9)

    ax_sau.set_xlabel("")
    ax_sau.set_ylabel("€/MWh", fontsize=10)
    ax_sau.grid(axis="x")
    ax_sau.set_ylim(0, 110)
    ax_sau.set_yticks(range(0, 120, 20))
    ax_sau.set_yticks(range(10, 120, 20), minor=True)
    for spine in ax_sau.spines.values():
        spine.set_visible(False)

    ax.annotate(
        "",
        xy=(0.655, 0.55),
        xytext=(0.62, 0.65),
        xycoords="axes fraction",
        arrowprops=dict(
            edgecolor=ARROW_COLOR,
            facecolor=ARROW_COLOR,
            linewidth=1.5,
            arrowstyle="-|>",
            connectionstyle="arc3,rad=0.2",
        ),
    )

    # inset: Australia steel imports

    ax_aus = ax.inset_axes([0.74, 0.15, 0.036, 0.25])

    composition_aus.T.plot.bar(ax=ax_aus, stacked=True, color=COLORS, legend=False)

    ax_aus.set_title("steel\nimports", fontsize=9)

    ax_aus.set_xlabel("")
    ax_aus.set_ylabel("€/tonne", fontsize=10)
    ax_aus.grid(axis="x")
    ax_aus.set_yticks(range(0, 601, 100))
    ax_aus.set_yticks(range(50, 601, 100), minor=True)
    ax_aus.set_ylim(0, 600)
    for spine in ax_aus.spines.values():
        spine.set_visible(False)

    ax.annotate(
        "",
        xy=(0.775, 0.35),
        xytext=(0.82, 0.315),
        xycoords="axes fraction",
        arrowprops=dict(
            edgecolor=ARROW_COLOR,
            facecolor=ARROW_COLOR,
            linewidth=1.5,
            arrowstyle="-|>",
            connectionstyle="arc3,rad=0.2",
        ),
    )

    # inset: land eligibility of Buenos Aires

    ax_land = ax.inset_axes([0.315, 0.08, 0.29, 0.29])

    shape.to_crs(crs.proj4_init).plot(
        ax=ax, color="none", edgecolor="k", linestyle=":", linewidth=1
    )

    add_land_eligibility_example(ax_land, shape, glc_fn, wdpa_fn)

    ax_land.set_title("wind exclusion\nzones (purple)", fontsize=9)

    ax.annotate(
        "",
        xy=(0.41, 0.22),
        xytext=(0.35, 0.17),
        xycoords="axes fraction",
        arrowprops=dict(
            edgecolor=ARROW_COLOR,
            facecolor=ARROW_COLOR,
            linewidth=1.5,
            arrowstyle="-|>",
            connectionstyle="arc3,rad=0.2",
        ),
    )

    for fn in snakemake.output:
        plt.savefig(fn, bbox_inches="tight")
