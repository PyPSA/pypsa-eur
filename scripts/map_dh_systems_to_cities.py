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
    """Mirror the label cleaning used in identify_district_heating_subnodes."""
    return (
        pd.Series([label]).fillna("DH").astype(str)
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
        usecols=["geonameid", "asciiname", "country_code", "latitude", "longitude", "population"],
    )
    df = df[df["country_code"].isin(countries)].copy()
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

    # Nearest city per system
    joined = gpd.sjoin_nearest(
        dh_areas.set_geometry("centroid"),
        cities,
        how="left",
        distance_col="distance_m",
    )

    joined["distance_km"] = joined["distance_m"] / 1000
    joined = joined.rename(columns={"asciiname": "city", "geonameid": "geonames_id"})

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
    joined[["latitude", "longitude"]] = joined.to_crs("EPSG:4326").geometry.apply(
        lambda p: pd.Series({"latitude": p.y, "longitude": p.x})
    )

    joined[out_cols].rename(columns={"Label": "label_original"}).to_csv(
        out_path, index=False
    )


if __name__ == "__main__":
    main(snakemake)  # type: ignore[name-defined]
