# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Maps the per-capita GDP and population values to non-NUTS3 regions.

The script takes as input the country code, a GeoDataFrame containing
the regions, and the file paths to the datasets containing the GDP and
POP values for non-NUTS3 countries.
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import rasterio
import xarray as xr
from _helpers import configure_logging, set_scenario_config
from rasterio.mask import mask
from shapely.geometry import box

logger = logging.getLogger(__name__)


def calc_gdp_pop(country, regions, gdp_non_nuts3, pop_non_nuts3):
    """
    Calculate the GDP p.c. and population values for non NUTS3 regions.

    Parameters:
    country (str): The two-letter country code of the non-NUTS3 region.
    regions (GeoDataFrame): A GeoDataFrame containing the regions.
    gdp_non_nuts3 (str): The file path to the dataset containing the GDP p.c values
    for non NUTS3 countries (e.g. MD, UA)
    pop_non_nuts3 (str): The file path to the dataset containing the POP values
    for non NUTS3 countries (e.g. MD, UA)

    Returns:
    tuple: A tuple containing two GeoDataFrames:
        - gdp: A GeoDataFrame with the mean GDP p.c. values mapped to each bus.
        - pop: A GeoDataFrame with the summed POP values mapped to each bus.
    """
    regions = regions.rename(columns={"name": "Bus"}).set_index("Bus")
    regions = regions[regions.country == country]
    # Create a bounding box for UA, MD from region shape, including a buffer of 10000 metres
    bounding_box = (
        gpd.GeoDataFrame(geometry=[box(*regions.total_bounds)], crs=regions.crs)
        .to_crs(epsg=3857)
        .buffer(10000)
        .to_crs(regions.crs)
    )

    # GDP Mapping
    logger.info(f"Mapping mean GDP p.c. to non-NUTS3 region: {country}")
    with xr.open_dataset(gdp_non_nuts3) as src_gdp:
        src_gdp = src_gdp.where(
            (src_gdp.longitude >= bounding_box.bounds.minx.min())
            & (src_gdp.longitude <= bounding_box.bounds.maxx.max())
            & (src_gdp.latitude >= bounding_box.bounds.miny.min())
            & (src_gdp.latitude <= bounding_box.bounds.maxy.max()),
            drop=True,
        )
        gdp = src_gdp.to_dataframe().reset_index()
    gdp = gdp.rename(columns={"GDP_per_capita_PPP": "gdp"})
    gdp = gdp[gdp.time == gdp.time.max()]
    gdp_raster = gpd.GeoDataFrame(
        gdp,
        geometry=gpd.points_from_xy(gdp.longitude, gdp.latitude),
        crs="EPSG:4326",
    )
    gdp_mapped = gpd.sjoin(gdp_raster, regions, predicate="within")
    gdp = (
        gdp_mapped.copy()
        .groupby(["Bus", "country"])
        .agg({"gdp": "mean"})
        .reset_index(level=["country"])
    )

    # Population Mapping
    logger.info(f"Mapping summed population to non-NUTS3 region: {country}")
    with rasterio.open(pop_non_nuts3) as src_pop:
        # Mask the raster with the bounding box
        out_image, out_transform = mask(src_pop, bounding_box, crop=True)
        out_meta = src_pop.meta.copy()
        out_meta.update(
            {
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
            }
        )
    masked_data = out_image[0]  # Use the first band (rest is empty)
    row_indices, col_indices = np.where(masked_data != src_pop.nodata)
    values = masked_data[row_indices, col_indices]

    # Affine transformation from pixel coordinates to geo coordinates
    x_coords, y_coords = rasterio.transform.xy(out_transform, row_indices, col_indices)
    pop_raster = pd.DataFrame({"x": x_coords, "y": y_coords, "pop": values})
    pop_raster = gpd.GeoDataFrame(
        pop_raster,
        geometry=gpd.points_from_xy(pop_raster.x, pop_raster.y),
        crs=src_pop.crs,
    )
    pop_mapped = gpd.sjoin(pop_raster, regions, predicate="within")
    pop = (
        pop_mapped.groupby(["Bus", "country"])
        .agg({"pop": "sum"})
        .reset_index()
        .set_index("Bus")
    )
    gdp_pop = regions.join(gdp.drop(columns="country"), on="Bus").join(
        pop.drop(columns="country"), on="Bus"
    )
    gdp_pop.fillna(0, inplace=True)

    return gdp_pop


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_gdp_pop_non_nuts3")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    regions = gpd.read_file(snakemake.input.regions)

    gdp_non_nuts3 = snakemake.input.gdp_non_nuts3
    pop_non_nuts3 = snakemake.input.pop_non_nuts3

    subset = {"MD", "UA"}.intersection(snakemake.params.countries)

    gdp_pop = pd.concat(
        [
            calc_gdp_pop(country, regions, gdp_non_nuts3, pop_non_nuts3)
            for country in subset
        ],
        axis=0,
    )

    logger.info(
        f"Exporting GDP and POP values for non-NUTS3 regions {snakemake.output}"
    )
    gdp_pop.reset_index().to_file(snakemake.output[0], driver="GeoJSON")
