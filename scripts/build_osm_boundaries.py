# SPDX-FileCopyrightText: : 2017-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import json
import logging

import country_converter as coco
import geopandas as gpd
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from build_shapes import eez
from shapely import line_merge
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon

logger = logging.getLogger(__name__)

cc = coco.CountryConverter()

GEO_CRS = "EPSG:4326"
EXCLUDER_LIST = [
    3788485,  # RU version of Sevastopol
    3795586,  # RU version of Crimea
]


def _create_linestring(row):
    """
    Create a LineString object from the given row.

    Parameters
    ----------
        row (dict): A dictionary containing the row data.

    Returns
    -------
        LineString: A LineString object representing the geometry.
    """
    coords = [(coord["lon"], coord["lat"]) for coord in row["geometry"]]
    return LineString(coords)


def _create_geometries(row, crs=GEO_CRS):
    """
    Create geometries from OSM data.

    This function processes OpenStreetMap (OSM) data to create geometries
    based on the roles of the members in the data. It supports "outer" and
    "inner" roles to construct polygons and multipolygons.

    Parameters
    ----------
        - row (dict): A dictionary containing OSM data with a "members" key.
        - crs (str, optional): Coordinate reference system for the geometries.
          Defaults to GEO_CRS.

    Returns
    -------
        - shapely.geometry.Polygon or shapely.geometry.MultiPolygon: The resulting
          geometry after processing the OSM data.
    """
    valid_roles = ["outer", "inner"]
    df = pd.json_normalize(row["members"])
    df = df[df["role"].isin(valid_roles) & ~df["geometry"].isna()]
    df.loc[:, "geometry"] = df.apply(_create_linestring, axis=1)

    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=crs)
    outer = line_merge(gdf[gdf["role"] == "outer"].union_all())

    if isinstance(outer, LineString):
        outer = Polygon(outer)
    if isinstance(outer, MultiLineString):
        polygons = [Polygon(line) for line in outer.geoms]
        outer = MultiPolygon(polygons)

    if not gdf[gdf["role"] == "inner"].empty:
        inner = line_merge(gdf[gdf["role"] == "inner"].union_all())

        if isinstance(inner, LineString):
            inner = Polygon(inner)
        if isinstance(inner, MultiLineString):
            polygons = [Polygon(line) for line in inner.geoms]
            inner = MultiPolygon(polygons)

        outer = outer.difference(inner)

    return outer


def build_osm_boundaries(country, adm1_path, offshore_shapes):
    """
    Build administrative boundaries from OSM data for a given country.

    Parameters
    ----------
        - country (str): The country code (e.g., 'DE' for Germany).
        - adm1_path (str): The file path to the administrative level 1 OSM data in JSON format.
        - offshore_shapes (GeoDataFrame): A GeoDataFrame containing offshore shapes to clip the boundaries.

    Returns
    -------
    GeoDataFrame: A GeoDataFrame containing the administrative boundaries with columns:
        - id: The administrative level 1 code.
        - country: The country code.
        - name: The administrative level 1 name.
        - osm_id: The OSM ID.
        - geometry: The geometry of the boundary.
    """
    logger.info(
        f"Building administrative boundaries from OSM data for country {country}."
    )

    df = json.load(open(adm1_path))
    df = pd.DataFrame(df["elements"])

    # Filter out ids in excluder list
    df = df[~df["id"].isin(EXCLUDER_LIST)]

    df.loc[:, "geometry"] = df.apply(_create_geometries, axis=1)

    col_tags = [
        "ISO3166-2",
        "name:en",
    ]

    tags = pd.json_normalize(df["tags"]).map(lambda x: str(x) if pd.notnull(x) else x)

    for ct in col_tags:
        if ct not in tags.columns:
            tags[ct] = pd.NA

    tags = tags.loc[:, col_tags]

    df = pd.concat([df, tags], axis="columns")
    df.drop(columns=["type", "tags", "bounds", "members"], inplace=True)

    # Rename columns
    df.rename(
        columns={
            "id": "osm_id",
            "ISO3166-2": "id",
            "name:en": "name",
        },
        inplace=True,
    )

    df["country"] = country

    # Resort columns
    df = df[["id", "country", "name", "osm_id", "geometry"]]

    # If any id is missing, sort by osm_id and number starting from 1
    if df["id"].isna().any():
        df.sort_values(by="osm_id", inplace=True)
        df["id"] = country + "-" + df.index.to_series().add(1).astype(str)

    gdf = gpd.GeoDataFrame(df, geometry="geometry", crs=GEO_CRS)

    # Clip gdf by offshore shapes
    gdf = gpd.overlay(gdf, offshore_shapes, how="difference")

    # Check if substring in "id" is equal to country, if not, drop
    gdf = gdf[gdf["id"].str.startswith(country)]

    return gdf


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_osm_boundaries", country="MD")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    country = snakemake.wildcards.country
    adm1_path = snakemake.input.json
    offshore_shapes = eez(snakemake.input.eez)

    boundaries = build_osm_boundaries(country, adm1_path, offshore_shapes)

    # Export
    boundaries.to_file(snakemake.output[0], driver="GeoJSON")
