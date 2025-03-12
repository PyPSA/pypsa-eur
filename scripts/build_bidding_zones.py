# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
from functools import reduce
from itertools import takewhile
from operator import attrgetter

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import MultiPolygon


def _simplify_polys(polys, minarea=0.1, tolerance=None, filterremote=True):
    if isinstance(polys, MultiPolygon):
        polys = sorted(polys.geoms, key=attrgetter("area"), reverse=True)
        mainpoly = polys[0]
        mainlength = np.sqrt(mainpoly.area / (2.0 * np.pi))
        if mainpoly.area > minarea:
            polys = MultiPolygon(
                [
                    p
                    for p in takewhile(lambda p: p.area > minarea, polys)
                    if not filterremote or (mainpoly.distance(p) < mainlength)
                ]
            )
        else:
            polys = mainpoly
    if tolerance is not None:
        polys = polys.simplify(tolerance=tolerance)
    return polys


def get_countries(naturalearth: str, country_list: set[str]) -> gpd.GeoDataFrame:
    """
    Load geometries for missing countries from Natural Earth dataset.

    Args:
        naturalearth: Path to Natural Earth dataset
        missing_countries: Set of country codes that are missing in the bidding zones

    Returns:
        GeoDataFrame: Contains geometries for all missing countries
    """
    df = gpd.read_file(naturalearth)

    # Names are a hassle in naturalearth, try several fields
    fieldnames = (
        df[x].where(lambda s: s != "-99") for x in ("ISO_A2", "WB_A2", "ADM0_A3")
    )
    df["name"] = reduce(lambda x, y: x.fillna(y), fieldnames, next(fieldnames)).str[:2]
    df.replace({"name": {"KV": "XK"}}, inplace=True)

    df = df.loc[
        df.name.isin(country_list) & ((df["scalerank"] == 0) | (df["scalerank"] == 5))
    ]
    return gpd.GeoDataFrame(
        df.set_index("name")[["geometry"]].geometry.map(_simplify_polys).set_crs(df.crs)
    )


def parse_zone_names(zone_names: pd.Series) -> tuple[set[str], pd.Series]:
    """
    Parse bidding zone names to extract country codes and handle combined zones.

    Args:
        zone_names: Series of zone names (e.g., 'DE', 'DE_LU', etc.)

    Returns:
        tuple containing:
            - set of country codes covered by the zones
            - Series of combined country strings (e.g., 'DE,LU')
    """
    # Extract primary country codes (first two characters)
    primary_countries = set(zone_names.str[:2])

    # Parse combined zones (e.g., 'DE_LU' -> 'LU')
    secondary_countries = (
        zone_names.str.split("_", expand=True)[1]
        .dropna()
        .loc[lambda x: x.str.len() == 2]
    )

    # Create combined country strings (e.g., 'DE,LU')
    country_strings = zone_names.str[:2] + (
        "," + secondary_countries.reindex(zone_names.index)
    ).fillna("")

    # Remove secondary countries from primary set
    covered_countries = primary_countries | set(secondary_countries)

    return covered_countries, country_strings


def split_combined_zones(gdf, country_shapes):
    """Split multi-country zones into individual country entries."""
    crs = gdf.crs
    new_rows = []

    for _, row in gdf.iterrows():
        countries = row.country.split(",")
        if len(countries) == 1:
            new_rows.append(row)
            continue

        # Split geometry for combined countries
        for country in countries:
            country_geom = country_shapes.loc[country].geometry
            intersection = row.geometry.intersection(country_geom)

            if not intersection.is_empty:
                new_row = row.copy()
                new_row.country = country
                new_row.geometry = intersection
                new_rows.append(new_row)

    return gpd.GeoDataFrame(new_rows, crs=crs)


def extract_shape_by_bbox(
    gdf: gpd.GeoDataFrame,
    country: str,
    min_lon: float,
    max_lon: float,
    min_lat: float,
    max_lat: float,
    region_id: str,
):
    """
    Extracts a shape from a country's GeoDataFrame based on latitude and longitude bounds.

    Parameters
    ----------
        - gdf (GeoDataFrame): GeoDataFrame containing country geometries.
        - country (str): The country code or name to filter.
        - min_lon, max_lon (float): Longitude bounds for extraction.
        - min_lat, max_lat (float): Latitude bounds for extraction.
        - region_id (str): String to assign an ID to the extracted region.

    Returns
    -------
        - gdf_new: Updated GeoDataFrame with the extracted shape separated.
    """
    country_gdf = gdf.explode().query(f"country == '{country}'").reset_index(drop=True)

    extracted_region = country_gdf.cx[min_lon:max_lon, min_lat:max_lat].assign(
        zone_name=region_id
    )

    remaining_country = (
        country_gdf.drop(extracted_region.index).dissolve(by="country").reset_index()
    )

    return pd.concat(
        [
            gdf.query(f"country != '{country}'"),
            remaining_country,
            extracted_region.dissolve(by="country").reset_index(),
        ]
    ).reset_index(drop=True)


def format_names(s: str):
    s = (
        s.replace("DK-DK1", "DKW1")
        .replace("DK-DK2", "DKE1")
        .replace("GB", "UK")
        .replace("IT_NORD", "ITN1")
        .replace("IT_SUD", "ITS1")
        .replace("LU", "LUG1")
        .replace("NO-NO3", "NOM1")
        .replace("NO-NO4", "NON1")
        .replace("SE-SE", "SE0")
        .replace("_", "")
        .ljust(4, "0")
    )[:4]
    return s


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_bidding_zones")

    # Load core bidding zones and country shapes
    countries = snakemake.params.countries

    rename_columns = {"zoneName": "zone_name", "countryKey": "country"}
    bidding_zones_elecmaps = gpd.read_file(
        snakemake.input.bidding_zones_electricitymaps
    )
    bidding_zones_elecmaps = bidding_zones_elecmaps.rename(columns=rename_columns)
    bidding_zones = bidding_zones_elecmaps[
        bidding_zones_elecmaps.country.isin(countries)
    ].drop(columns="countryName")

    if not set(countries).issubset(bidding_zones.country):
        raise ValueError("Missing countries in electricitymaps bidding zones")

    # Manual corrections: replace Italy by entsoepy version
    bidding_zones = bidding_zones[bidding_zones.country != "IT"]

    bidding_zones_entsoe = gpd.read_file(snakemake.input.bidding_zones_entsoepy)
    bidding_zones_entsoe = bidding_zones_entsoe.rename(
        columns={"zoneName": "zone_name"}
    )
    core_countries, country_strings = parse_zone_names(
        bidding_zones_entsoe["zone_name"]
    )
    bidding_zones_entsoe["country"] = country_strings
    italian_zones = bidding_zones_entsoe[bidding_zones_entsoe.country == "IT"]
    bidding_zones = pd.concat([bidding_zones, italian_zones], ignore_index=True)
    bidding_zones = bidding_zones.sort_values("zone_name").reset_index(drop=True)

    # manual corrections: remove islands
    islands = [
        "DK-BHM",
        "ES-CE",
        "ES-CN-HI",
        "ES-CN-IG",
        "ES-CN-LP",
        "ES-CN-LZ",
        "ES-IB-FO",
        "ES-IB-IZ",
        "ES-IB-ME",
        "ES-IB-MA",
        "ES-CN-FV",
        "ES-CN-GC",
        "ES-CN-TE",
        "ES-ML",
        "GB-ORK",
        "GB-ZET",
        "PT-MA",
        "PT-AC",
    ]
    bidding_zones = bidding_zones[~bidding_zones.zone_name.isin(islands)]

    # manually merge southern norwegian zones
    nos0_idx = bidding_zones.query("zone_name in ['NO-NO1', 'NO-NO2', 'NO-NO5']").index
    bidding_zones = pd.concat(
        [
            bidding_zones.drop(nos0_idx),
            bidding_zones.loc[nos0_idx]
            .dissolve(by="country")
            .reset_index()
            .assign(zone_name="NOS0"),
        ]
    )

    # Extract Crete
    bidding_zones = extract_shape_by_bbox(
        bidding_zones,
        country="GR",
        min_lon=24.0,
        max_lon=26.5,
        min_lat=35.0,
        max_lat=35.7,
        region_id="GR03",
    )

    # manually add zone group for DE_LU
    bidding_zones.loc[
        bidding_zones.zone_name.isin(["DE", "LU"]), "cross_country_zone"
    ] = "DE_LU"

    # if turkey is not in the list of countries, drop northern cyprus
    if "TR" not in countries:
        bidding_zones = bidding_zones[~bidding_zones.zone_name.eq("XX")]

    # rename zones
    bidding_zones["zone_name"] = bidding_zones["zone_name"].apply(format_names)

    bidding_zones.to_file(snakemake.output.file)
