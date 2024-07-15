# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""

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


def calc_gdp_ppp(country, regions, gdp_non_nuts3, ppp_non_nuts3):
    """
    Calculate the GDP and PPP values for non NUTS3 regions.

    Parameters:
    country (str): The two-letter country code of the non-NUTS3 region.
    regions (GeoDataFrame): A GeoDataFrame containing the regions.
    gdp_non_nuts3 (str): The file path to the dataset containing the GDP values
    for non NUTS3 countries (e.g. MD, UA)
    ppp_non_nuts3 (str): The file path to the dataset containing the PPP values
    for non NUTS3 countries (e.g. MD, UA)

    Returns:
    tuple: A tuple containing two GeoDataFrames:
        - gdp: A GeoDataFrame with the aggregated GDP values mapped to each bus.
        - ppp: A GeoDataFrame with the aggregated PPP values mapped to each bus.
    """
    regions = regions.drop(columns=["x", "y"])
    regions = regions[regions.country == country]
    # Create a bounding box for UA, MD from region shape, including a buffer of 10000 metres
    bounding_box = (
        gpd.GeoDataFrame(geometry=[box(*regions.total_bounds)], crs=regions.crs)
        .to_crs(epsg=3857)
        .buffer(10000)
        .to_crs(regions.crs)
    )

    # GDP
    logger.info(f"Mapping GDP values to non-NUTS3 region: {regions.country.unique()}")
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
    gdp = gpd.GeoDataFrame(
        gdp,
        geometry=gpd.points_from_xy(gdp.longitude, gdp.latitude),
        crs="EPSG:4326",
    )
    gdp = gpd.sjoin(gdp, regions, predicate="within")
    gdp = (
        gdp.groupby(["Bus", "country"])
        .agg({"gdp": "sum"})
        .reset_index(level=["country"])
    )

    # PPP
    logger.info(f"Mapping PPP values to non-NUTS3 region: {regions.country.unique()}")
    with rasterio.open(ppp_non_nuts3) as src_ppp:
        # Mask the raster with the bounding box
        out_image, out_transform = mask(src_ppp, bounding_box, crop=True)
        out_image,
        out_meta = src_ppp.meta.copy()
        out_meta.update(
            {
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
            }
        )
    masked_data = out_image[0]  # Use the first band (rest is empty)
    row_indices, col_indices = np.where(masked_data != src_ppp.nodata)
    values = masked_data[row_indices, col_indices]

    # Affine transformation from pixel coordinates to geo coordinates
    x_coords, y_coords = rasterio.transform.xy(out_transform, row_indices, col_indices)
    ppp = pd.DataFrame({"x": x_coords, "y": y_coords, "ppp": values})
    ppp = gpd.GeoDataFrame(
        ppp,
        geometry=gpd.points_from_xy(ppp.x, ppp.y),
        crs=src_ppp.crs,
    )
    ppp = gpd.sjoin(ppp, regions, predicate="within")
    ppp = (
        ppp.groupby(["Bus", "country"])
        .agg({"ppp": "sum"})
        .reset_index()
        .set_index("Bus")
    )
    gdp_ppp = regions.join(gdp.drop(columns="country"), on="Bus").join(
        ppp.drop(columns="country"), on="Bus"
    )
    gdp_ppp.fillna(0, inplace=True)

    return gdp_ppp


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_gdp_ppp_non_nuts3")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    n = pypsa.Network(snakemake.input.base_network)
    substation_lv_i = n.buses.index[n.buses["substation_lv"]]
    regions = (
        gpd.read_file(snakemake.input.regions)
        .set_index("name")
        .reindex(substation_lv_i)
    )

    gdp_non_nuts3 = snakemake.input.gdp_non_nuts3
    ppp_non_nuts3 = snakemake.input.ppp_non_nuts3

    countries_non_nuts3 = pd.Index(("MD", "UA"))
    subset = set(countries_non_nuts3) & set(snakemake.params.countries)

    gdp_ppp = pd.concat(
        [
            calc_gdp_ppp(country, regions, gdp_non_nuts3, ppp_non_nuts3)
            for country in subset
        ],
        axis=0,
    )

    logger.info(
        f"Exporting GDP and PPP values for non-NUTS3 regions {snakemake.output}"
    )
    gdp_ppp.reset_index().to_file(snakemake.output, driver="GeoJSON")
