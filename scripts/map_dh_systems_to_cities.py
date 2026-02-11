# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create a static lookup of district-heating systems to nearby city names.

Inputs
------
- dh_areas (GeoPackage/GeoJSON): polygons with at least columns
  * country (ISO2)
  * Label (system name/id)
- geonames_cities (TSV): downloaded Geonames cities500 file

Output
------
CSV with columns:
  label_original, subnode_label, country, city, geonames_id,
  population, distance_km, latitude, longitude

Notes
-----
- No online geocoder calls at runtime; Geonames dataset is downloaded once
  in the retrieve stage.
- Distances are centroid-to-city-point; CRS projected to EPSG:3857.
"""

import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd

logger = logging.getLogger(__name__)


def sanitize_label(label: str) -> str:
    """Mirror the label cleaning used in build_district_heating_subnodes."""
    return (
        pd.Series([label])
        .fillna("DH")
        .astype(str)
        .str.replace(r"[^\w\s-]", "", regex=True)
        .str.strip()
        .str.replace(r"\s+", "_", regex=True)
    ).iloc[0]


def load_geonames(path: str, countries: set[str]) -> gpd.GeoDataFrame:
    cols = [
        "geonameid",
        "name",
        "asciiname",
        "alternatenames",
        "latitude",
        "longitude",
        "feature_class",
        "feature_code",
        "country_code",
        "cc2",
        "admin1",
        "admin2",
        "admin3",
        "admin4",
        "population",
        "elevation",
        "dem",
        "timezone",
        "modification_date",
    ]
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        usecols=[
            "geonameid",
            "asciiname",
            "country_code",
            "latitude",
            "longitude",
            "population",
            "feature_code",
        ],
    )
    # Filter to requested countries
    df = df[df["country_code"].isin(countries)].copy()

    # Filter to actual cities only (PPL* but not PPLX which are districts/neighborhoods)
    # PPL = populated place, PPLA-PPLA4 = admin capitals at various levels
    # PPLX = section of populated place (district), PPLL = populated locality
    city_codes = ["PPL", "PPLA", "PPLA2", "PPLA3", "PPLA4", "PPLA5", "PPLC"]
    df = df[df["feature_code"].isin(city_codes)].copy()

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["longitude"], df["latitude"]),
        crs="EPSG:4326",
    )
    return gdf


def main(snakemake):
    dh_path = Path(snakemake.input.dh_areas)
    cities_path = Path(snakemake.input.cities)
    out_path = Path(snakemake.output[0])
    out_path.parent.mkdir(parents=True, exist_ok=True)

    dh_areas = gpd.read_file(dh_path)
    if dh_areas.empty:
        pd.DataFrame(
            columns=[
                "label_original",
                "subnode_label",
                "country",
                "city",
                "geonames_id",
                "population",
                "distance_km",
                "latitude",
                "longitude",
            ]
        ).to_csv(out_path, index=False)
        return

    if "country" not in dh_areas.columns:
        raise ValueError("dh_areas must contain a 'country' column")
    if "Label" not in dh_areas.columns:
        raise ValueError("dh_areas must contain a 'Label' column for system id")

    dh_areas = dh_areas.to_crs("EPSG:3857")
    dh_areas["centroid"] = dh_areas.geometry.centroid
    dh_areas["subnode_label"] = dh_areas["Label"].apply(sanitize_label)

    countries = set(dh_areas["country"].unique())
    cities = load_geonames(cities_path, countries).to_crs("EPSG:3857")

    if cities.empty:
        logger.warning("No cities found for countries %s", countries)
        cities["distance"] = pd.NA
        cities["asciiname"] = pd.NA

    # Strategy: Find the LARGEST city associated with each DH area.
    # 1. Check cities within a buffered DH area polygon (10 km buffer in EPSG:3857)
    # 2. Among candidates, pick the one with the highest population
    # 3. Fall back to nearest city if no candidates found
    #
    # The buffer ensures we capture the primary city even when its centre-point
    # lies just outside the DH network polygon boundary (common for suburban
    # DH systems like Leverkusen vs Cologne).
    BUFFER_M = 10_000  # 10 km buffer in metres (EPSG:3857)
    results = []

    for idx, dh_row in dh_areas.iterrows():
        dh_geom = dh_row.geometry
        dh_centroid = dh_row["centroid"]

        # Find cities within the buffered DH area polygon
        buffered_geom = dh_geom.buffer(BUFFER_M)
        cities_nearby = cities[cities.geometry.within(buffered_geom)]

        if not cities_nearby.empty:
            # Select the largest city by population within the buffered polygon
            largest_city = cities_nearby.loc[cities_nearby["population"].idxmax()]
            distance_m = dh_centroid.distance(largest_city.geometry)
            city_name = largest_city["asciiname"]
            geonames_id = largest_city["geonameid"]
            population = largest_city["population"]
            lat = largest_city["latitude"]
            lon = largest_city["longitude"]
        else:
            # Fall back to nearest city
            distances = cities.geometry.distance(dh_centroid)
            if distances.empty:
                city_name = None
                geonames_id = None
                population = None
                distance_m = None
                lat = None
                lon = None
            else:
                nearest_idx = distances.idxmin()
                nearest_city = cities.loc[nearest_idx]
                distance_m = distances.loc[nearest_idx]
                city_name = nearest_city["asciiname"]
                geonames_id = nearest_city["geonameid"]
                population = nearest_city["population"]
                lat = nearest_city["latitude"]
                lon = nearest_city["longitude"]

        results.append(
            {
                "Label": dh_row["Label"],
                "subnode_label": dh_row["subnode_label"],
                "country": dh_row["country"],
                "city": city_name,
                "geonames_id": geonames_id,
                "population": population,
                "distance_km": distance_m / 1000 if distance_m is not None else None,
                "latitude": lat,
                "longitude": lon,
            }
        )

    joined = pd.DataFrame(results)

    out_cols = [
        "Label",
        "subnode_label",
        "country",
        "city",
        "geonames_id",
        "population",
        "distance_km",
        "latitude",
        "longitude",
    ]

    joined[out_cols].rename(columns={"Label": "label_original"}).to_csv(
        out_path, index=False
    )


if __name__ == "__main__":
    main(snakemake)  # type: ignore[name-defined]
