# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build the capacities of existing cement, steel, ammonia and methanol plants per model region in t/a.

Operating integrated cement plants come from the Global Cement and Concrete Tracker. Operating blast furnaces and direct reduction furnaces come from the Global Iron and Steel Tracker, with the most recent relining as build year and the announced retirement as phase-out year. Operating ammonia plants and methanol plants with fossil feedstock come from the Global Chemicals Tracker, with the start of the plant as build year. Each plant is assigned to the nearest model region; plants further than 100 km away from any region are dropped. Plants with unknown capacity get the 20% quantile of the capacities of their carrier.

References
----------
- [Global Energy Monitor: Global Iron and Steel Tracker](https://globalenergymonitor.org/projects/global-iron-and-steel-tracker/)
- [Global Energy Monitor: Global Chemicals Tracker](https://globalenergymonitor.org/projects/global-chemicals-tracker/)
- [Global Energy Monitor: Global Cement and Concrete Tracker](https://globalenergymonitor.org/projects/global-cement-and-concrete-tracker/)
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pycountry

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)


def country_to_code(name):
    try:
        return pycountry.countries.lookup(name).alpha_2
    except LookupError:
        return None


def map_to_regions(
    gdf: gpd.GeoDataFrame, regions: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Assign plants to the nearest model region and drop plants further than 100 km away."""
    crs = "EPSG:3035"
    gdf = gdf.to_crs(crs).sjoin_nearest(
        regions.to_crs(crs), how="inner", max_distance=100e3
    )
    return gdf.rename(columns={"name": "bus"})


def prepare_gem_database(regions):
    """
    Load GEM database of cement plants and map onto bus regions.
    """
    df = pd.read_excel(
        snakemake.input.gem_gcct,
        sheet_name="Final data",
        na_values=["N/A", "unknown", ">0"],
    )

    df["country_code"] = df["Country/area"].apply(country_to_code)
    # XK not in pycountry
    df.loc[df["Country/area"] == "Kosovo", "country_code"] = "XK"

    df = df[df.country_code.isin(snakemake.params.countries)]

    df["Start date"] = pd.to_numeric(
        df["Start date"].astype(str).str[:4], errors="coerce"
    )

    latlon = (
        df["Coordinates"]
        .str.split(", ", expand=True)
        .rename(columns={0: "lat", 1: "lon"})
    )
    geometry = gpd.points_from_xy(latlon["lon"], latlon["lat"])
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

    gdf = map_to_regions(gdf, regions)
    gdf["country"] = gdf.bus.str[:2]

    # operating plants with a kiln; grinding plants only process clinker
    status = ["operating", "operating pre-retirement"]
    gdf = gdf[gdf["Operating status"].isin(status)]
    gdf = gdf[gdf["Plant type"] == "integrated"]

    gdf.rename(
        columns={
            "Cement capacity (million metric tonnes per annum)": "p_set",
            "Start date": "build_year",
        },
        inplace=True,
    )
    gdf["p_set"] *= 1e6
    gdf["carrier"] = "cement"
    gdf["Out"] = 0

    cement = gdf[["bus", "country", "carrier", "p_set", "build_year", "Out"]]

    return cement


def prepare_gem_iron_units(
    regions: gpd.GeoDataFrame,
    iron_units: str,
    steel_plants: str,
) -> pd.DataFrame:
    """
    Load GEM database of blast furnaces and direct reduction furnaces and map onto bus regions.

    The build year is the most recent relining or, if unknown, the start of
    the unit. The phase-out year is the announced retirement.

    Parameters
    ----------
    regions : gpd.GeoDataFrame
        onshore regions on which the units are mapped to
    iron_units : str
        path to the iron unit data of the GEM Global Iron and Steel Tracker
    steel_plants : str
        path to the plant-level data of the tracker with the coordinates
    """
    kwargs = dict(na_values=["N/A", "unknown", ">0"])
    sheets = {"BOF": "Blast furnaces", "gas DRI": "Direct reduced iron furnaces"}
    df = pd.concat(
        {
            carrier: pd.read_excel(iron_units, sheet_name=sheet, **kwargs)
            for carrier, sheet in sheets.items()
        },
        names=["carrier"],
    ).reset_index("carrier")
    df = df[df["Unit status"].isin(["operating", "operating pre-retirement"])]

    df["country"] = df["Country/area"].apply(country_to_code)
    # XK not in pycountry
    df.loc[df["Country/area"] == "Kosovo", "country"] = "XK"
    df = df[df.country.isin(snakemake.params.countries)]

    coordinates = pd.read_excel(steel_plants, sheet_name="Plant data", index_col=0)
    latlon = (
        coordinates.loc[df["GEM plant ID"], "Coordinates"]
        .str.split(", ", expand=True)
        .astype(float)
    )
    geometry = gpd.points_from_xy(latlon[1], latlon[0])
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")
    gdf = map_to_regions(gdf, regions)

    def year(dates: pd.Series) -> pd.Series:
        return pd.to_numeric(dates.astype(str).str[:4], errors="coerce")

    gdf["build_year"] = year(gdf["Most recent relining"]).fillna(
        year(gdf["Start date"])
    )
    gdf["Out"] = year(gdf["Retired date"])
    # kt/a to t/a
    gdf["p_set"] = gdf["Current capacity (ttpa)"] * 1e3

    return gdf[["bus", "country", "carrier", "p_set", "build_year", "Out"]]


def prepare_gem_chemicals(regions: gpd.GeoDataFrame, chemicals: str) -> pd.DataFrame:
    """
    Load GEM database of ammonia and fossil methanol plants and map onto bus regions.

    Parameters
    ----------
    regions : gpd.GeoDataFrame
        onshore regions on which the plants are mapped to
    chemicals : str
        path to the GEM Global Chemicals Tracker
    """
    plants = pd.read_excel(
        chemicals, sheet_name="Plant Data", index_col=0, na_values=["not found"]
    )
    df = pd.read_excel(
        chemicals, sheet_name="Chemical Capacities", na_values=["not found"]
    ).join(plants, on="GEM Plant ID")
    df = df[df["Capacity Status"].isin(["operating", "operating pre-retirement"])]

    # methanol plants without fossil feedstock are not grey methanol plants
    fossil = df["Feedstock"].str.contains("natural gas|crude oil|naphtha|coal")
    carriers = {"ammonia": "Haber-Bosch", "methanol": "grey methanol"}
    df["carrier"] = df["Product"].map(carriers)
    df = df[df.carrier.notna() & (fossil | (df.Product == "ammonia"))]

    df["country"] = df["Country/Area"].apply(country_to_code)
    # XK not in pycountry
    df.loc[df["Country/Area"] == "Kosovo", "country"] = "XK"
    df = df[df.country.isin(snakemake.params.countries)]

    latlon = df["Coordinates"].str.split(", ", expand=True).astype(float)
    geometry = gpd.points_from_xy(latlon[1], latlon[0])
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")
    gdf = map_to_regions(gdf, regions)

    gdf["build_year"] = gdf["Start Date | Year"]
    gdf["Out"] = np.nan
    # kt/a to t/a
    gdf["p_set"] = gdf["Capacity (ttpa)"] * 1e3

    return gdf[["bus", "country", "carrier", "p_set", "build_year", "Out"]]


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_industry_plants")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    # get cement plant data
    cement = prepare_gem_database(regions)
    # steel plants
    steel = prepare_gem_iron_units(
        regions, snakemake.input.gem_iron_units, snakemake.input.gem_gspt
    )
    # ammonia and methanol plants
    chemicals = prepare_gem_chemicals(regions, snakemake.input.gem_gct)

    industry_plants = pd.concat([cement, steel, chemicals], ignore_index=True)

    # conservative capacity for plants with unknown capacity
    fill = industry_plants.groupby("carrier").p_set.transform("quantile", 0.2)
    industry_plants["p_set"] = industry_plants.p_set.fillna(fill)

    industry_plants.to_csv(snakemake.output.industry_plants)
