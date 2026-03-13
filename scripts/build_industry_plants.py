# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Adds existing power and heat generation capacities for initial planning
horizon.
"""

import logging
from types import SimpleNamespace

import country_converter as coco
import geopandas as gpd
import numpy as np
import pandas as pd
import pycountry

from scripts._helpers import (
    configure_logging,
    set_scenario_config,
)

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()
idx = pd.IndexSlice
spatial = SimpleNamespace()


def country_to_code(name):
    try:
        return pycountry.countries.lookup(name).alpha_2
    except LookupError:
        return None


def prepare_gem_database(regions):
    """
    Load GEM database of cement plants and map onto bus regions.
    """
    df = pd.read_excel(
        snakemake.input.gem_gcct,
        sheet_name="Plant Data",
        na_values=["N/A", "unknown", ">0"],
    )

    df["country_code"] = df["Country/Area"].apply(country_to_code)

    df = df[df.country_code.isin(snakemake.params.countries)]

    df["Start date"] = pd.to_numeric(
        df["Start date"].str.split("-").str[0], errors="coerce"
    )

    latlon = (
        df["Coordinates"]
        .str.split(", ", expand=True)
        .rename(columns={0: "lat", 1: "lon"})
    )
    geometry = gpd.points_from_xy(latlon["lon"], latlon["lat"])
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

    gdf = gpd.sjoin(gdf, regions, how="inner", predicate="within")

    gdf.rename(columns={"name": "bus"}, inplace=True)
    gdf["country"] = gdf.bus.str[:2]

    # just get bus, country, grouping_year, plant_size, date_in, date_out,
    # get plants that are still operating
    gdf = gdf[gdf["Operating status"] == "operating"]

    gdf.rename(
        columns={
            "Cement Capacity (millions metric tonnes per annum)": "p_set",
            "Start date": "build_year",
        },
        inplace=True,
    )
    gdf["p_set"] *= 1e6
    gdf["carrier"] = "cement"
    gdf["Out"] = 0

    cement = gdf[["bus", "country", "carrier", "p_set", "build_year", "Out"]]

    return cement


def prepare_plant_data(
    regions: gpd,
    isi_database: str,
) -> tuple[pd.DataFrame, gpd.GeoDataFrame]:
    """
    Reads in the Fraunhofer ISI database with high resolution plant data and maps them to the bus regions.
    Returns the database as df.

    Parameters
    ----------
    regions : gpd
        onshore regions on which the industry sites are mapped to
    isi_database: str
        path to the fraunhofer isi database
    """

    isi_data = pd.read_excel(isi_database, sheet_name="Database", index_col=1)
    # assign bus region to each plant
    geometry = gpd.points_from_xy(isi_data["Longitude"], isi_data["Latitude"])
    plant_data = gpd.GeoDataFrame(isi_data, geometry=geometry, crs="EPSG:4326")
    plant_data = gpd.sjoin(plant_data, regions, how="inner", predicate="within")
    plant_data.rename(columns={"name": "bus"}, inplace=True)
    # filter for countries in model scope
    plant_data = plant_data[plant_data.Country.isin(snakemake.params.countries)]
    # replace UK with GB in Country column
    plant_data["country"] = plant_data["Country"].replace("UK", "GB")
    # add carrier column
    carrier_dict = {
        "Blast furnace": "BOF",
        "Direct reduction NG": "gas DRI",
        "Ammonia SMR": "Haber-Bosch",
        "Methanol SMR": "grey methanol",
    }
    plant_data["carrier"] = plant_data["Process status qup"].replace(carrier_dict)
    # get build_year
    plant_data["build_year"] = (
        plant_data["Year of last modernisation"]
        .replace("x", np.nan)
        .fillna(plant_data["Last Relining"])
    )
    plant_data.rename(
        columns={"Production in tons (calibrated)": "p_set"}, inplace=True
    )
    plant_data = plant_data[["bus", "country", "carrier", "p_set", "build_year", "Out"]]

    return plant_data


def prepare_ammonia_data(regions, plant_data):
    df = pd.read_csv(snakemake.input.ammonia, index_col=0)

    geometry = gpd.points_from_xy(df.Longitude, df.Latitude)
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

    gdf = gpd.sjoin(gdf, regions, how="inner", predicate="within")

    gdf.rename(columns={"name": "bus"}, inplace=True)
    gdf["country"] = gdf.bus.str[:2]
    # filter for countries that are missing
    gdf = gdf[
        (~gdf.country.isin(plant_data.country.unique()))
        & (gdf.country.isin(snakemake.params.countries))
    ]
    # following approach from build_industrial_distribution_key.py
    for country in gdf.Country:
        facilities = gdf.query("country == @country")
        production = facilities["Ammonia [kt/a]"]
        # assume 50% of the minimum production for missing values
        production = production.fillna(0.5 * facilities["Ammonia [kt/a]"].min())

    # missing data
    gdf.drop(gdf[gdf["Ammonia [kt/a]"].isna()].index, inplace=True)

    # get average plant age:
    avg_age = plant_data[plant_data.carrier == "Haber-Bosch"]["build_year"].mean()
    # restructure to fit other data sources
    gdf["build_year"] = avg_age
    gdf["Out"] = 0
    gdf["carrier"] = "Haber-Bosch"
    gdf["Ammonia [kt/a]"] *= 1e3
    gdf.rename(columns={"Ammonia [kt/a]": "p_set"}, inplace=True)
    gdf = gdf[["bus", "country", "carrier", "p_set", "build_year", "Out"]]

    return gdf


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industry_plants",
            configfiles="config/config.default.yaml",
            clusters="39",
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    regions = gpd.read_file(snakemake.input.regions_onshore).set_index("name")

    # get cement plant data
    cement = prepare_gem_database(regions)
    # get data from fraunhofer database
    fh_data = prepare_plant_data(regions, snakemake.input.isi_database)
    # ammonia plants non-EU27
    ammonia = prepare_ammonia_data(regions, fh_data)

    industry_plants = pd.concat([cement, fh_data, ammonia], ignore_index=True)

    industry_plants.to_csv(snakemake.output.industry_plants)
