# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Creates GIS shape files of NUTS3 and OSM ADM1 areas (for BA, MD, UA, and XK),
clipped to onshore territory and enriched with GDP and population data.
Optionally assigns bidding zones when administrative clustering is enabled.
"""

import logging
import unicodedata

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.mask import mask
from shapely.geometry import box

from scripts._helpers import _simplify_polys, configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


GEO_CRS = "EPSG:4326"
DISTANCE_CRS = "EPSG:3035"
GDP_YEAR = 2019
POP_YEAR = 2019
DROP_REGIONS = [
    "ES703",
    "ES704",
    "ES705",
    "ES706",
    "ES707",
    "ES708",
    "ES709",
    "ES630",
    "ES640",
    "FRY10",
    "FRY20",
    "FRY30",
    "FRY40",
    "FRY50",
    "NO0B1",
    "NO0B2",
    "PT200",
    "PT300",
]
OTHER_GDP_TOTAL_2019 = {  # in bn. USD
    "BA": 20.48,  # World Bank
    "MD": 11.74,  # World Bank
    "UA": 153.9,  # World Bank
    "XK": 7.9,  # https://de.statista.com/statistik/daten/studie/415738/umfrage/bruttoinlandsprodukt-bip-des-kosovo/
}
OTHER_POP_2019 = {  # in 1000 persons
    "BA": 3361,  # World Bank
    "MD": 2664,  # World Bank
    "UA": 44470,  # World Bank
    "XK": 1782,  # World Bank
}
EXCHANGE_EUR_USD_2019 = 1.1
NUTS3_INCLUDE = [
    "DE80N",
    "DEF08",
    "DK014",
    "DK050",
    "GBH34",
    "GBJ34",
    "GBJ43",
    "NL33A",
    "NL33C",
    "NL342",
    "NL411",
    "NL412",
    "PL428",
]


def normalise_text(text):
    text = unicodedata.normalize("NFD", text)
    text = "".join(char for char in text if unicodedata.category(char) != "Mn")
    text = text.replace("*", "")
    return "".join(char for char in text if char.isascii())


def simplify_europe(regions):
    logger.info(
        "Simplifying geometries for Europe by removing small islands smaller than 500 km2 or further than 200 km away."
    )
    coverage = (
        regions.to_crs(DISTANCE_CRS)
        .groupby("country")["geometry"]
        .apply(lambda x: x.union_all())
    )
    coverage_dk = coverage.loc[["DK"]]
    coverage = coverage.apply(_simplify_polys, minarea=500 * 1e6, maxdistance=200 * 1e3)
    coverage_dk = coverage_dk.apply(
        _simplify_polys, minarea=65 * 1e6, maxdistance=200 * 1e3
    )
    coverage.loc["DK"] = coverage_dk.values[0]
    coverage = gpd.GeoDataFrame(geometry=coverage, crs=DISTANCE_CRS).to_crs(GEO_CRS)

    coverage = pd.concat([coverage, regions.loc[NUTS3_INCLUDE, ["geometry"]]])
    shape = coverage.union_all()

    regions_polygon = regions.explode()
    regions_polygon = gpd.sjoin(
        regions_polygon,
        gpd.GeoDataFrame(geometry=[shape], crs=GEO_CRS),
        how="inner",
        predicate="intersects",
    )
    regions_polygon = regions_polygon.groupby(["level3"])["geometry"].apply(
        lambda x: x.union_all()
    )

    regions = regions.loc[regions.index.isin(regions_polygon.index)]
    regions.loc[regions_polygon.index, "geometry"] = regions_polygon
    return regions


def calc_gdp_pop(country, regions, gdp_non_nuts3, pop_non_nuts3):
    region = regions.loc[regions.country == country, ["geometry"]]
    bounding_box = (
        gpd.GeoDataFrame(geometry=[box(*region.total_bounds)], crs=region.crs)
        .to_crs(epsg=3857)
        .buffer(10000)
        .to_crs(region.crs)
    )

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
    gdp_mapped = gpd.sjoin(gdp_raster, region, predicate="within")
    gdp = gdp_mapped.groupby(["id"]).agg({"gdp": "mean"})

    logger.info(f"Mapping summed population to non-NUTS3 region: {country}")
    with rasterio.open(pop_non_nuts3) as src_pop:
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
    masked_data = out_image[0]
    row_indices, col_indices = np.where(masked_data != src_pop.nodata)
    values = masked_data[row_indices, col_indices]

    x_coords, y_coords = rasterio.transform.xy(out_transform, row_indices, col_indices)
    pop_raster = gpd.GeoDataFrame(
        {"x": x_coords, "y": y_coords, "pop": values},
        geometry=gpd.points_from_xy(x_coords, y_coords),
        crs=src_pop.crs,
    )
    pop_mapped = gpd.sjoin(pop_raster, region, predicate="within")
    pop = pop_mapped.groupby(["id"]).agg({"pop": "sum"}).reset_index().set_index("id")
    gdp_pop = region.join(gdp).join(pop).drop(columns="geometry")
    gdp_pop.fillna(0, inplace=True)

    gdp_pop["gdp"] = gdp_pop["gdp"].round(0)
    gdp_pop["pop"] = gdp_pop["pop"].div(1e3).round(0)
    gdp_pop["pop"] = (
        gdp_pop["pop"].div(gdp_pop["pop"].sum()).mul(OTHER_POP_2019[country]).round(0)
    )
    gdp_pop["gdp"] = (
        gdp_pop["gdp"]
        .mul(1e9)
        .div(gdp_pop["gdp"].sum())
        .mul(OTHER_GDP_TOTAL_2019[country])
        .div(EXCHANGE_EUR_USD_2019)
        / (1e3 * gdp_pop["pop"])
    ).round(0)
    return gdp_pop


def bidding_zone_map(
    regions: gpd.GeoDataFrame, bidding_zones: gpd.GeoDataFrame
) -> pd.Series:
    """
    Map bidding zones to regions on a country-by-country basis, assigning each region
    to the bidding zone with which it has the largest overlap. If a region doesn't
    overlap with any bidding zone, it's assigned to the nearest one as a fallback.

    Parameters
    ----------
    regions : geopandas.GeoDataFrame
        The regions GeoDataFrame with a 'country' column
    bidding_zones : geopandas.GeoDataFrame
        The bidding zones GeoDataFrame with a 'zone_name' and "country" column.

    Returns
    -------
    bidding_zone_map : pandas.Series
        A Series with the same index as regions, containing the assigned bidding zone names.
    """
    bz_map = pd.Series(pd.NA, regions.index, dtype="object")
    logger.info("Mapping bidding zones to regions based on maximum overlap area")

    for country in regions["country"].unique():
        country_regions = regions.loc[regions["country"] == country].copy()
        country_bz = bidding_zones[bidding_zones.country == country]

        if country_regions.empty or country_bz.empty:
            logger.debug(f"Skipping country {country}: no regions or bidding zones found")
            continue

        if country_regions.crs != country_bz.crs:
            country_bz = country_bz.to_crs(country_regions.crs)

        assignments, no_overlap = [], []
        for idx, region in country_regions.iterrows():
            best_zone, max_overlap = None, 0
            for _, zone in country_bz.iterrows():
                intersection = region.geometry.intersection(zone.geometry)
                if not intersection.is_empty and intersection.area > max_overlap:
                    max_overlap = intersection.area
                    best_zone = zone.zone_name
            if best_zone is not None:
                assignments.append((idx, best_zone))
            else:
                no_overlap.append(idx)

        for idx, zone in assignments:
            bz_map[idx] = zone

        if no_overlap:
            logger.info(
                f"Assigning {len(no_overlap)} regions in {country} to nearest bidding zone"
            )
            for idx in no_overlap:
                region_geom = regions.loc[idx, "geometry"]
                distances = country_bz.geometry.distance(region_geom)
                if not distances.empty:
                    bz_map[idx] = country_bz.loc[distances.idxmin(), "zone_name"]
                else:
                    logger.warning(f"Could not find any bidding zone for region {idx} in {country}")

    unassigned = bz_map.isnull().sum()
    if unassigned > 0:
        logger.warning(f"{unassigned} regions couldn't be assigned to any bidding zone")
    return bz_map


def create_regions(
    country_list,
    nuts3_path,
    ba_adm1_path,
    md_adm1_path,
    ua_adm1_path,
    xk_adm1_path,
    offshore_shapes,
    nuts3_gdp,
    nuts3_pop,
    bidding_zones_path,
    other_gdp,
    other_pop,
):
    """
    Build the NUTS3/ADM1 region GeoDataFrame enriched with GDP, population,
    and optionally bidding zone assignments.

    Parameters
    ----------
    country_list : list[str]
    nuts3_path : str
    ba_adm1_path, md_adm1_path, ua_adm1_path, xk_adm1_path : str
    offshore_shapes : geopandas.GeoDataFrame
        Used to clip onshore region boundaries at the coastline.
    nuts3_gdp, nuts3_pop : str
    bidding_zones_path : str or list
    other_gdp, other_pop : str

    Returns
    -------
    geopandas.GeoDataFrame
    """
    logger.info("Processing NUTS regions.")
    regions = gpd.read_file(nuts3_path)
    regions.loc[regions.CNTR_CODE == "EL", "CNTR_CODE"] = "GR"
    regions["NUTS_ID"] = regions["NUTS_ID"].str.replace("EL", "GR")
    regions.loc[regions.CNTR_CODE == "UK", "CNTR_CODE"] = "GB"
    regions["NUTS_ID"] = regions["NUTS_ID"].str.replace("UK", "GB")

    regions = regions[["NUTS_ID", "CNTR_CODE", "NAME_LATN", "geometry"]].rename(
        columns={"NUTS_ID": "id", "CNTR_CODE": "country", "NAME_LATN": "name"}
    )
    regions["id"] = regions["id"].apply(normalise_text)
    regions["level1"] = regions["id"].str[:3]
    regions["level2"] = regions["id"].str[:4]
    regions["level3"] = regions["id"]

    logger.info("Processing non-NUTS regions.")
    regions_non_nuts = pd.concat(
        [gpd.read_file(p) for p in [ba_adm1_path, md_adm1_path, ua_adm1_path, xk_adm1_path]]
    ).drop(columns=["osm_id"])
    regions_non_nuts["id"] = regions_non_nuts["id"].apply(normalise_text)
    regions_non_nuts["name"] = regions_non_nuts["name"].apply(normalise_text)
    regions_non_nuts["level1"] = regions_non_nuts["id"]
    regions_non_nuts["level2"] = regions_non_nuts["id"]
    regions_non_nuts["level3"] = regions_non_nuts["id"]

    regions["geometry"] = regions["geometry"].difference(
        regions_non_nuts.geometry.union_all()
    )

    logger.info("Harmonising NUTS and non-NUTS regions.")
    regions = pd.concat([regions, regions_non_nuts]).set_index("id")
    regions = regions.drop(DROP_REGIONS, errors="ignore")

    logger.info("Clipping regions by offshore shapes.")
    regions["geometry"] = regions["geometry"].difference(
        offshore_shapes.geometry.union_all()
    )

    logger.info(f"Importing JRC ARDECO GDP data for year {GDP_YEAR}.")
    nuts3_gdp_df = pd.read_csv(nuts3_gdp, index_col=[0])
    nuts3_gdp_df = nuts3_gdp_df.query("LEVEL_ID == 3 and UNIT == 'EUR'")
    nuts3_gdp_df.index = nuts3_gdp_df.index.str.replace("UK", "GB").str.replace("EL", "GR")
    regions["gdp"] = nuts3_gdp_df[str(GDP_YEAR)]

    logger.info(f"Importing JRC ARDECO population data for year {POP_YEAR}.")
    nuts3_pop_df = pd.read_csv(nuts3_pop, index_col=[0])
    nuts3_pop_df = nuts3_pop_df.query("LEVEL_ID == 3")
    nuts3_pop_df.index = nuts3_pop_df.index.str.replace("UK", "GB").str.replace("EL", "GR")
    regions["pop"] = nuts3_pop_df[str(POP_YEAR)].div(1e3).round(0)

    other_countries = {"BA", "MD", "UA", "XK"}
    if any(c in country_list for c in other_countries):
        gdp_pop = pd.concat(
            [calc_gdp_pop(c, regions, other_gdp, other_pop) for c in other_countries],
            axis=0,
        )
        regions.loc[gdp_pop.index, ["gdp", "pop"]] = gdp_pop[["gdp", "pop"]]

    regions = regions[
        ["name", "level1", "level2", "level3", "gdp", "pop", "country", "geometry"]
    ]
    regions.index.name = "index"
    regions = simplify_europe(regions)
    regions = regions.query("country in @country_list")

    if bidding_zones_path:
        bidding_zones = gpd.read_file(bidding_zones_path)
        regions = regions.assign(bidding_zone=bidding_zone_map(regions, bidding_zones))

    return regions


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_nuts3_shapes")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).set_index("name")

    regions = create_regions(
        snakemake.params.countries,
        snakemake.input.nuts3_2021,
        snakemake.input.ba_adm1,
        snakemake.input.md_adm1,
        snakemake.input.ua_adm1,
        snakemake.input.xk_adm1,
        offshore_shapes,
        snakemake.input.nuts3_gdp,
        snakemake.input.nuts3_pop,
        snakemake.input.bidding_zones,
        snakemake.input.other_gdp,
        snakemake.input.other_pop,
    )

    logger.info(
        f"Exporting NUTS3 and ADM1 shapes with GDP and POP values to {snakemake.output.nuts3_shapes}."
    )
    regions.reset_index().to_file(snakemake.output.nuts3_shapes)
