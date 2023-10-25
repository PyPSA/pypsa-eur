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
from rasterio.features import geometry_mask
from rasterio.plot import show

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()

AREA_CRS = "ESRI:54009"

NICE_NAMES = {
    "pipeline-h2": r"H$_2$ (pipeline)",
    "shipping-lh2": "H$_2$ (ship)",
    "shipping-ftfuel": "Fischer-Tropsch",
    "shipping-meoh": "methanol",
    "shipping-lch4": "methane",
    "shipping-lnh3": "ammonia",
    "shipping-steel": "steel",
}

def rename(s):
    if "solar" in s:
        return "solar"
    if "wind" in s:
        return "wind"
    if "storage" in s or "inverter" in s:
        return "storage"
    if "transport" in s or "shipping fuel"  in s or "dry bulk" in s or "pipeline" in s:
        return "transport"
    if "evaporation" in s or "liquefaction" in s or "compress" in s:
        return "evaporation/liquefaction"
    if "direct air capture" in s or "heat pump" in s:
        return "direct air capture"
    if s in ["Haber-Bosch (exp)", "Fischer-Tropsch (exp)", "methanolisation (exp)", "methanation (exp)"]:
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
    composition = df.query(query_str).groupby(["esc", "subcategory", "importer"]).value.min()

    minimal = {}
    for name, group in composition.groupby("esc"):
        c = group.unstack("importer").droplevel("esc")
        minimal[name] = c[c.sum().idxmin()]
    composition = pd.concat(minimal, axis=1)

    composition = composition.groupby(rename).sum().div(production)

    composition = composition.where(composition > 0.01).dropna(how='all')

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
            layer="WDPA_Oct2023_Public_shp-polygons",
        ).to_crs(AREA_CRS)
        excluder.add_geometry(wdpa.geometry, buffer=1000)

    band, transform = shape_availability(shape, excluder)
    mask = ~geometry_mask(
        [shape.geometry.values[0]],
        transform=transform,
        invert=False,
        out_shape=band.shape
    )
    masked_band = np.where(mask, ~band, np.nan)

    shape.plot(ax=ax, color="none", edgecolor='k', linewidth=1)
    show(masked_band, transform=transform, cmap="Purples", ax=ax)
    ax.set_axis_off()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_import_world_map",
            simpl="",
            opts="",
            clusters="100",
            ll="v1.5",
            sector_opts="Co2L0-2190SEG-T-H-B-I-S-A-imp",
            planning_horizons="2050",
            configfiles="../../config/config.100n-seg.yaml",
        )

    configure_logging(snakemake)

    import_fn = snakemake.input.imports
    profile_fn = snakemake.input.profiles
    gadm_fn = snakemake.input.gadm_arg[0]
    glc_fn = snakemake.input.copernicus_glc[0]
    wdpa_fn = "/home/fneum/bwss/playgrounds/pr/pypsa-eur/resources/WDPA_Oct2023.gpkg"

    config = snakemake.config

    plt.style.use(['bmh', snakemake.input.rc])

    tech_colors = config["plotting"]["tech_colors"]

    COLORS = {
        "wind": tech_colors["onwind"],
        "solar": tech_colors["solar"],
        "storage": tech_colors["battery"],
        "electrolysis": tech_colors["H2 Electrolysis"],
        "direct air capture": tech_colors["DAC"],
        "hydrogen conversion": tech_colors["Fischer-Tropsch"],
        "iron ore": "#4e4f55",
        "direct iron reduction": tech_colors["steel"],
        "electric arc furnace": "#8795a8",
        "evaporation/liquefaction": "#8487e8",
        "transport": "#f7a572",
    }

    # load capacity factor time series

    ds = xr.open_dataset(profile_fn)
    profile = ds.sel(exporter='MA', importer='EUE').p_max_pu.to_pandas().T
    profile.rename(columns={"onwind": "wind", "solar-utility": "solar"}, inplace=True)

    # download world country shapes, version clipped to Europe, and GADM in AR

    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres')).set_index("iso_a3")
    world.drop("ATA", inplace=True)

    eu_countries = cc.convert(config["countries"], src="iso2", to="iso3")
    europe = world.loc[eu_countries]

    gadm = gpd.read_file(gadm_fn, layer="ADM_ADM_1").set_index("NAME_1")
    shape = gadm.to_crs(AREA_CRS).loc[["Buenos Aires"]].geometry

    # load import costs

    df = pd.read_csv(import_fn, sep=";", keep_default_na=False)

    # bugfix for Namibia
    df["exporter"] = df.exporter.replace("", "NA")

    import_costs = df.query("subcategory == 'Cost per MWh delivered' and esc == 'shipping-meoh'").groupby("exporter").value.min()
    import_costs.index = cc.convert(import_costs.index.str.split("-").str[0], src='iso2', to='iso3')

    import_costs.drop("RUS", inplace=True, errors="ignore")

    composition_arg = get_cost_composition(
        df,
        "AR",
        ["shipping-lh2", "shipping-ftfuel", "shipping-meoh", "shipping-lch4", "shipping-lnh3"],
        500e6
    )

    composition_sau = get_cost_composition(
        df,
        "SA",
        ["pipeline-h2", "shipping-lh2"],
        500e6
    )

    composition_aus = get_cost_composition(
        df,
        "AU",
        ["shipping-steel"],
        100e6
    )

    # add to legend

    composition_arg.loc["iron ore"] = pd.NA
    composition_arg.loc["direct iron reduction"] = pd.NA
    composition_arg.loc["electric arc furnace"] = pd.NA

    # create plot

    crs = ccrs.EqualEarth()

    fig, ax = plt.subplots(figsize=(14,14), subplot_kw={"projection": crs})

    # main axis: choropleth layer

    world.to_crs(crs).plot(
        column=import_costs.reindex(world.index),
        linewidth=1,
        edgecolor="black",
        ax=ax,
        cmap='Greens_r',
        legend=True,
        vmin=100,
        vmax=150,
        legend_kwds=dict(label="Cost for methanol fuel delivered to Europe [€/MWh]", orientation="horizontal", extend='max', shrink=.6, aspect=30, pad=.01),
        missing_kwds=dict(color="#eee", label="not considered"),
    )

    europe.to_crs(crs).plot(
        linewidth=1,
        edgecolor="black",
        ax=ax,
        color='lavender',
    )

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_facecolor('none')
    fig.set_facecolor('none')

    ax.text(0.93, 0.01, "Projection:\nEqual Earth", transform=ax.transAxes, fontsize=9, color="grey")

    plt.tight_layout()

    # inset: wind and solar profiles

    ax_prof = ax.inset_axes([0.36, 0.68, 0.09, 0.11])

    week_profile = profile.loc["2013-03-01":"2013-03-07", ["solar", "wind"]]

    week_profile.plot(
        ax=ax_prof,
        linewidth=1,
        color=['gold', "royalblue"],
        ylim=(0,1),
        clip_on=False,
    )

    ax_prof.legend(
        title="",
        loc=(0,1),
        fontsize=8,
        ncol=2,
        columnspacing=0.8
    )
    ax_prof.set_xlabel("March 2013", fontsize=8)
    ax_prof.set_ylabel("profile [p.u.]", fontsize=8)
    ax_prof.tick_params(axis='both', labelsize=8)

    ax_prof.xaxis.set_major_locator(mdates.DayLocator())
    ax_prof.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    xticks = week_profile.resample('D').mean().index
    ax_prof.set_xticks(xticks)
    ax_prof.set_xticklabels(xticks.day)

    for spine in ax_prof.spines.values():
        spine.set_visible(False)

    for label in ax_prof.get_xticklabels():
        label.set_fontsize(8)

    ax.annotate(
        '', 
        xy=(0.45, 0.75),
        xytext=(0.485, 0.72),
        xycoords='axes fraction',
        arrowprops=dict(
        edgecolor='#555',
        facecolor='#555',
        linewidth=1.5,
        arrowstyle='-|>',
        connectionstyle="arc3,rad=0.2"
        )
    )

    # inset: Argentina e-fuel import costs

    ax_arg = ax.inset_axes([0.07, 0.08, 0.09, 0.33])

    composition_arg.T.plot.bar(ax=ax_arg, stacked=True, color=COLORS)

    handles, labels = ax_arg.get_legend_handles_labels()
    handles.reverse()
    labels.reverse()

    ax_arg.legend(handles, labels, title="", ncol=1, fontsize=9, loc=(1,0))

    ax_arg.set_title("Import costs from\nArgentina to Europe", fontsize=9)

    ax_arg.set_xlabel("")
    ax_arg.set_ylim(0, 110)
    ax_arg.set_yticks(range(0, 111, 20))
    ax_arg.set_yticks(range(10, 111, 20), minor=True)
    ax_arg.set_ylabel("€/MWh", fontsize=10)
    ax_arg.grid(axis="x")
    for spine in ax_arg.spines.values():
        spine.set_visible(False)

    ax.annotate(
        '', 
        xy=(0.25, 0.15),
        xytext=(0.33, 0.2),
        xycoords='axes fraction',
        arrowprops=dict(
            edgecolor='#555',
            facecolor='#555',
            linewidth=1.5,
            arrowstyle='-|>',
            connectionstyle="arc3,rad=-0.2"
        )
    )

    # inset: Saudi Arabia hydrogen pipeline versus ship imports

    ax_sau = ax.inset_axes([0.6625, 0.33, 0.036, 0.2])

    composition_sau.T.plot.bar(ax=ax_sau, stacked=True, color=COLORS, legend=False)

    ax_sau.set_title(r"LH$_2$" + " ship\nvs. pipeline", fontsize=9)

    ax_sau.set_xlabel("")
    ax_sau.set_ylabel("€/MWh", fontsize=10)
    ax_sau.grid(axis="x")
    ax_sau.set_ylim(0, 90)
    ax_sau.set_yticks(range(0, 91, 20))
    ax_sau.set_yticks(range(10, 91, 20), minor=True)
    for spine in ax_sau.spines.values():
        spine.set_visible(False)

    ax.annotate(
        '', 
        xy=(0.655, 0.55),
        xytext=(0.62, 0.65),
        xycoords='axes fraction',
        arrowprops=dict(
            edgecolor='#555',
            facecolor='#555',
            linewidth=1.5,
            arrowstyle='-|>',
            connectionstyle="arc3,rad=0.2"
        )
    )

    # inset: Australia steel imports

    ax_aus = ax.inset_axes([0.75, 0.15, 0.018, 0.25])

    composition_aus.T.plot.bar(ax=ax_aus, stacked=True, color=COLORS, legend=False)

    ax_aus.set_title("steel\nimports", fontsize=9)

    ax_aus.set_xlabel("")
    ax_aus.set_ylabel("€/tonne", fontsize=10)
    ax_aus.grid(axis="x")
    ax_aus.set_yticks(range(0, 600, 100))
    ax_aus.set_yticks(range(50, 600, 100), minor=True)
    for spine in ax_aus.spines.values():
        spine.set_visible(False)


    ax.annotate(
        '', 
        xy=(0.77, 0.35),
        xytext=(0.815, 0.31),
        xycoords='axes fraction',
        arrowprops=dict(
            edgecolor='#555',
            facecolor='#555',
            linewidth=1.5,
            arrowstyle='-|>',
            connectionstyle="arc3,rad=0.2"
        )
    )

    # inset: land eligibility of Buenos Aires

    ax_land = ax.inset_axes([0.315, 0.08, 0.29, 0.29])

    shape.to_crs(crs.proj4_init).plot(ax=ax, color="none", edgecolor='k', linestyle=":", linewidth=1)

    add_land_eligibility_example(ax_land, shape, glc_fn, wdpa_fn)

    ax_land.set_title("wind exclusion\nzones (purple)", fontsize=9)

    ax.annotate(
        '', 
        xy=(0.41, 0.22),
        xytext=(0.35, 0.17),
        xycoords='axes fraction',
        arrowprops=dict(
            edgecolor='#555',
            facecolor='#555',
            linewidth=1.5,
            arrowstyle='-|>',
            connectionstyle="arc3,rad=0.2"
        )
    )

    for fn in snakemake.output:
        plt.savefig(fn, bbox_inches="tight")
