# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Combines bidding zone shape files from two sources. The `electricitymaps-contrib` data is more accurate and are used as the baseline. The Italian bidding zones from `entsoe-py` are more preferred and are used to override the baseline. Manual adjustments are made to match the TYNDP 2024 configuration. Small islands are removed and Crete is considered as independent of Greece. Southern Norwegian zones are merged.

Outputs
-------

- ``resources/bidding_zones.geojson``:
"""

import geopandas as gpd
import pandas as pd


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

    return country_strings


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
    country_gdf = gdf.explode().query("country == @country").reset_index(drop=True)

    extracted_region = country_gdf.cx[min_lon:max_lon, min_lat:max_lat].assign(
        zone_name=region_id
    )

    remaining_country = (
        country_gdf.drop(extracted_region.index).dissolve(by="country").reset_index()
    )

    return pd.concat(
        [
            gdf.query("country != @country"),
            remaining_country,
            extracted_region.dissolve(by="country").reset_index(),
        ]
    ).reset_index(drop=True)


def format_names(s: str):
    s = (
        s.replace("DK-DK1", "DKW1")
        .replace("DK-DK2", "DKE1")
        .replace("FR-C", "FR15")
        .replace("GB", "UK")
        .replace("UK-N", "UKNI")
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
    country_strings = parse_zone_names(bidding_zones_entsoe["zone_name"])
    bidding_zones_entsoe["country"] = country_strings
    italian_zones = bidding_zones_entsoe[bidding_zones_entsoe.country == "IT"]
    bidding_zones = pd.concat([bidding_zones, italian_zones], ignore_index=True)

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

    bidding_zones = bidding_zones.sort_values("zone_name").reset_index(drop=True)
    bidding_zones.to_file(snakemake.output.file)
