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
from shapely.geometry import MultiPolygon, Polygon


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


def replace_country(
    source: gpd.GeoDataFrame,
    reference: gpd.GeoDataFrame,
    country: str,
    default_tolerance: float = 0.05,
    tolerance_dict: dict[str, dict[str, float]] = None,
):
    """
    Replace the shape of a specified country in the source shapes file with the corresponding shape from a reference shapes file.

    Parameters
    ----------
    source : geopandas.GeoDataFrame
        Original shapes file, including the country to replace.
    reference : geopandas.GeoDataFrame
        Alternative shapes file to use for the replaced the country.
    country : str
        The country code to be replaced in the source shapes file.
    default_tolerance : float, optional
        Default snapping tolerance (in degrees) used to align neighboring borders.
    tolerance_dict : dict of dict, optional
        A nested dictionary specifying custom tolerances (in degrees) for snapping between specific zone pairs.
        Format: {zone_name: {neighbor_zone_name: tolerance}}.
        If not provided, the default tolerance value is used.

    Returns
    -------
    geopandas.GeoDataFrame
        A shapes file with the specified country replaced and geometries snapped to
        neighbors to avoid overlaps.
    """

    # Remove country from the source
    bidding_zones = source[source.country != country]

    # Get the data from reference and add to source
    country_strings = parse_zone_names(reference["zone_name"])
    reference["country"] = country_strings
    country_zones = reference[reference.country == country]
    bidding_zones = pd.concat([bidding_zones, country_zones], ignore_index=True)

    # Loop on zones in source country
    for z in bidding_zones.query("country==@country").zone_name:
        # Find neighbors
        z_geom = bidding_zones.loc[bidding_zones.zone_name == z].iloc[0].geometry
        neighbors = bidding_zones[
            (bidding_zones.intersects(z_geom)) & (bidding_zones.country != country)
        ]

        if neighbors.empty:
            continue

        # Snap borders for each neighbor
        for n in neighbors.zone_name:
            zi = bidding_zones.query("zone_name == @z")
            ni = bidding_zones.query("zone_name == @n")
            tol = tolerance_dict.get(z, default_tolerance).get(n, default_tolerance)
            ni.loc[:, "geometry"] = (
                ni.snap(zi, tolerance=tol, align=False)
                .buffer(0)
                # the new neighbor border cannot overlap the reference zone
                .difference(zi, align=False)
            )
            bidding_zones = pd.concat(
                [bidding_zones.query("zone_name != @n"), ni], ignore_index=True
            )

        # Remove neighbors overlaps
        for n in neighbors.zone_name:
            ni = bidding_zones.query("zone_name == @n")
            ni.loc[:, "geometry"] = ni.difference(
                neighbors.query("zone_name!=@n").dissolve(), align=False
            )
            bidding_zones = pd.concat(
                [bidding_zones.query("zone_name != @n"), ni], ignore_index=True
            )

    return bidding_zones


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


def remove_holes(geom):
    if geom.geom_type == "Polygon":
        return Polygon(geom.exterior)
    elif geom.geom_type == "MultiPolygon":
        return MultiPolygon([Polygon(p.exterior) for p in geom.geoms])
    else:
        return geom


def format_names(s: str):
    s = (
        s.replace("DK-DK1", "DKW1")
        .replace("DK-DK2", "DKE1")
        .replace("ES-CN", "ES")
        .replace("ES-IB", "ES")
        .replace("FR-C", "FR15")
        .replace("UK-N", "UKNI")
        .replace("UK", "GB")
        .replace("IT_NORD", "ITN1")
        .replace("IT_SUD", "ITS1")
        .replace("LU", "LUG1")
        .replace("NO-NO1", "NOS1")
        .replace("NO-NO2", "NOS2")
        .replace("NO-NO3", "NOM1")
        .replace("NO-NO4", "NON1")
        .replace("NO-NO5", "NOS5")
        .replace("SE-SE", "SE0")
        .replace("_", "")
        .replace("-", "")
        .ljust(4, "0")
    )[:4]
    return s


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_bidding_zones", configfiles="config/test/config.clusters.yaml"
        )

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

    bidding_zones_entsoe = gpd.read_file(snakemake.input.bidding_zones_entsoepy)
    bidding_zones_entsoe = bidding_zones_entsoe.rename(
        columns={"zoneName": "zone_name"}
    )

    tolerance_dict = {
        "IT_NORD": {
            "FR": 0.05,
            "CH": 0.04,
            "AT": 0.024,
            "SI": 0.01,
        }
    }

    if "IT" in countries:
        bidding_zones = replace_country(
            source=bidding_zones,
            reference=bidding_zones_entsoe,
            country="IT",
            tolerance_dict=tolerance_dict,
        )

    if snakemake.params.remove_islands:
        # manual corrections: remove islands
        islands = [
            # Bornholm
            "DK-BHM",
            # Canary Islands
            "ES-CN-HI",
            "ES-CN-IG",
            "ES-CN-LP",
            "ES-CN-LZ",
            "ES-CN-FV",
            "ES-CN-GC",
            "ES-CN-TE",
            # Balearic Islands
            "ES-IB-FO",
            "ES-IB-IZ",
            "ES-IB-ME",
            "ES-IB-MA",
            # Melilla & Ceuta
            "ES-ML",
            "ES-CE",
            # Orkney Islands
            "GB-ORK",
            # Shetland Islands
            "GB-ZET",
            # Madeira & Azores Islands
            "PT-MA",
            "PT-AC",
        ]
        bidding_zones = bidding_zones[~bidding_zones.zone_name.isin(islands)]

    if snakemake.params.aggregate_to_tyndp:
        # Manually merge southern norwegian zones
        nos0_idx = bidding_zones.query(
            "zone_name in ['NO-NO1', 'NO-NO2', 'NO-NO5']"
        ).index
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

    # remove holes from geometries
    bidding_zones["geometry"] = bidding_zones["geometry"].apply(remove_holes)

    # if turkey is not in the list of countries, drop northern cyprus
    if "TR" not in countries:
        bidding_zones = bidding_zones[~bidding_zones.zone_name.eq("XX")]

    # rename zones
    bidding_zones["zone_name"] = bidding_zones["zone_name"].apply(format_names)

    bidding_zones = bidding_zones.sort_values("zone_name").reset_index(drop=True)
    bidding_zones.to_file(snakemake.output.file)
