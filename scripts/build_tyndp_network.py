# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import logging

import geopandas as gpd
import pandas as pd
from _helpers import configure_logging, set_scenario_config
from shapely.geometry import LineString, Point

logger = logging.getLogger(__name__)

GEO_CRS = "EPSG:4326"
DISTANCE_CRS = "EPSG:3035"
BUSES_COLUMNS = [
    "station_id",
    "voltage",
    "dc",
    "symbol",
    "under_construction",
    "tags",
    "x",
    "y",
    "country",
    "geometry",
]
LINES_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "circuits",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
LINKS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "length",
    "underground",
    "under_construction",
    "tags",
    "geometry",
]
TRANSFORMERS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage_bus0",
    "voltage_bus1",
    "s_nom",
    "station_id",
    "geometry",
]
CONVERTERS_COLUMNS = [
    "bus0",
    "bus1",
    "voltage",
    "p_nom",
    "geometry",
]


def format_bz_names(s: str):
    s = s.replace("FR-C", "FR15").replace("UK-N", "UKNI").replace("UK", "GB")
    return s


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
        id=region_id
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


def build_shapes(bz_fn, geo_crs: str = GEO_CRS, distance_crs: str = DISTANCE_CRS):
    """
    Process bidding zones from the shape file and calculate representative point. Deduce the country shapes and their representative point.

    Parameters
    ----------
        - bz_fn (str | Path): Path to bidding zone shape file.
        - geo_crs (CRS, optional): Coordinate reference system for geographic calculations. Defaults to GEO_CRS.
        - distance_crs (CRS, optional): Coordinate reference system to use for distance calculations. Defaults to DISTANCE_CRS.

    Returns
    -------
        - bidding_shapes: A GeoDataFrame including bidding zone geometry, representative point and id.
        - country_shapes: A GeoDataFrame including country geometry and representative point.
    """
    bidding_zones = gpd.read_file(bz_fn)

    # Bidding zone shapes
    bidding_shapes = bidding_zones.assign(
        bz_id=lambda df: df["zone_name"].apply(format_bz_names),
        node=lambda df: df.geometry.to_crs(distance_crs)
        .representative_point()
        .to_crs(geo_crs),
        x=lambda df: df["node"].x,
        y=lambda df: df["node"].y,
    ).set_index("bz_id")

    # Country shapes
    country_shapes = bidding_shapes.dissolve(by="country")[["geometry"]].assign(
        node=lambda df: df.geometry.to_crs(distance_crs)
        .representative_point()
        .to_crs(geo_crs),
        x=lambda df: df["node"].x,
        y=lambda df: df["node"].y,
    )

    # Correct DK, IT, GR and SE coordinates
    country_shapes.loc["DK", ["node", "x", "y"]] = bidding_shapes.loc[
        "DKE1", ["node", "x", "y"]
    ]
    country_shapes.loc["IT", ["node", "x", "y"]] = bidding_shapes.loc[
        "ITCA", ["node", "x", "y"]
    ]
    country_shapes.loc["GR", ["node", "x", "y"]] = bidding_shapes.loc[
        "GR00", ["node", "x", "y"]
    ]
    country_shapes.loc["SE", ["node", "x", "y"]] = bidding_shapes.loc[
        "SE01", ["node", "x", "y"]
    ]

    return bidding_shapes, country_shapes


def build_buses(
    buses_fn,
    bidding_shapes: gpd.GeoDataFrame,
    country_shapes: gpd.GeoDataFrame,
    geo_crs: str = GEO_CRS,
):
    """
    Extend the node list for both electricity and hydrogen with attributes, incl. country and coordinates.

    Parameters
    ----------
        - buses_fn (str | Path): Path to bidding zone shape file.
        - bidding_shapes (GeoDataFrame): A GeoDataFrame including bidding zone geometry, representative point and id.
        - country_shapes (GeoDataFrame): A GeoDataFrame including country geometry and representative point.
        - geo_crs (CRS, optional): Coordinate reference system for geographic calculations. Defaults to GEO_CRS.


    Returns
    -------
        - buses: A GeoDataFrame of electrical buses including country and coordinates.
        - buses_h2: A GeoDataFrame of hydrogen buses including country and coordinates.
    """
    buses = (
        pd.read_excel(buses_fn)
        .merge(
            bidding_shapes[["country", "node", "x", "y"]],
            how="outer",
            left_on="NODE",
            right_index=True,
        )
        .rename({"NODE": "bus_id", "node": "geometry"}, axis=1)
        .assign(
            station_id=lambda df: df["bus_id"],
            voltage=380,  # TODO Improve assumption
            dc=None,
            symbol="Substation",
            under_construction="f",
            tags=lambda df: df["bus_id"],
        )
        .set_index("bus_id")[BUSES_COLUMNS]
    )
    buses = gpd.GeoDataFrame(buses, geometry="geometry", crs=geo_crs)

    # Assume the same coordinates for all LU buses
    buses.loc["LUB1"] = buses.loc["LUB1"].fillna(buses.loc["LUG1"])
    buses.loc["LUF1"] = buses.loc["LUF1"].fillna(buses.loc["LUG1"])
    buses.loc["LUV1"] = buses.loc["LUV1"].fillna(buses.loc["LUG1"])

    # Manually add Italian virtual nodes  # TODO Refine assumptions
    buses.loc["ITCO"] = (
        buses.loc[["FR15"]]
        .assign(station_id="ITCO", country="IT", tags="ITCO")
        .loc["FR15"]
    )
    buses.loc["ITVI"] = (
        buses.loc[["ITSI"]].assign(station_id="ITVI", tags="ITVI").loc["ITSI"]
    )

    buses_h2 = (
        country_shapes[["node", "x", "y"]]
        .reset_index()
        .rename({"node": "geometry"}, axis=1)
        .assign(
            bus_id=lambda df: df[["country"]] + " H2",
            station_id=lambda df: df["bus_id"],
            voltage=None,
            dc="f",
            symbol="Substation",
            under_construction="f",
            tags=lambda df: df["bus_id"],
        )
        .set_index("bus_id")[BUSES_COLUMNS]
    )
    buses_h2 = gpd.GeoDataFrame(buses_h2, geometry="geometry", crs=geo_crs)

    # Manually add IBIT and IBFI nodes  # TODO Refine assumptions
    buses_h2.loc["IBIT H2"] = (
        buses.loc[["ITN1"]]
        .assign(station_id="IBIT H2", voltage=None, dc="f", tags="IBIT H2")
        .loc["ITN1"]
    )
    ibfi_lat, ibfi_long = 63.0, 25.0
    buses_h2.loc["IBFI H2"] = (
        buses_h2.loc[["FI H2"]]
        .assign(
            station_id="IBFI H2",
            tags="IBFI H2",
            x=ibfi_long,
            y=ibfi_lat,
            geometry=Point(ibfi_long, ibfi_lat),
        )
        .loc["FI H2"]
    )

    return buses, buses_h2


def format_grid_names(s: str):
    s = (
        s
        # Poland organizes its lines in three sections,
        # PL00 for demand/generation, -E for exporting lines and -I for importing lines
        .replace("PL00E", "PL00")
        .replace("PL00I", "PL00")
        .replace("UK", "GB")
    )
    return s


def build_links(
    grid_fn,
    buses: gpd.GeoDataFrame,
    geo_crs: str = GEO_CRS,
    distance_crs: str = DISTANCE_CRS,
):
    """
    Process reference grid information to produce link data. p_nom are NTC values.

    Parameters
    ----------
        - grid_fn (str | Path): Path to bidding zone shape file.
        - geo_crs (CRS, optional): Coordinate reference system for geographic calculations. Defaults to GEO_CRS.
        - distance_crs (CRS, optional): Coordinate reference system to use for distance calculations. Defaults to DISTANCE_CRS.

    Returns
    -------
        - links: A GeoDataFrame including NTC from the reference grid.

    """
    links = pd.read_excel(grid_fn)
    links["Border"] = links["Border"].apply(format_grid_names)
    links[["bus0", "bus1"]] = links.Border.str.split("-", expand=True)

    # Create forward and reverse direction dataframes
    # TODO: combine to bidirectional links
    forward_links = links[["bus0", "bus1", "Summary Direction 1"]].rename(
        columns={"Summary Direction 1": "p_nom"}
    )

    reverse_links = links[["bus1", "bus0", "Summary Direction 2"]].rename(
        columns={"bus1": "bus0", "bus0": "bus1", "Summary Direction 2": "p_nom"}
    )

    # Combine into unidirectional links
    links = pd.concat([forward_links, reverse_links])

    # Add missing attributes
    links = links.merge(
        buses["geometry"], how="left", left_on="bus0", right_index=True
    ).merge(
        buses["geometry"],
        how="left",
        left_on="bus1",
        right_index=True,
        suffixes=("0", "1"),
    )

    unknown_buses = set(
        links["bus0"][links[["bus0", "geometry0"]].isna().any(axis=1)]
    ).union(set(links["bus1"][links[["bus1", "geometry1"]].isna().any(axis=1)]))
    known_exceptions = {
        "DEKF",  # Connection from DE to the Kriegers Flak offshore wind farm
        "DKKF",  # Connection from DK to the Kriegers Flak offshore wind farm
        "DZ00",  # Algeria
        "EG00",  # Egypt
        "IS00",  # Iceland
        "IL00",  # Israel
        "LY00",  # Libya
        "MA00",  # Morocco
        "MD00",  # Moldova
        "PS00",  # Palestine
        "TN00",  # Tunisia
        "TR00",  # Turkey
        "UA00",  # Ukraine
        "UA01",  # Ukraine
    }
    if unknown_buses - known_exceptions:
        logger.warning(
            f"Dropping links connected to unknown buses: "
            f"{', '.join(sorted(unknown_buses - known_exceptions))}"
        )
    links = links.dropna()  # TODO Remove this when all nodes are known

    links["geometry"] = gpd.GeoSeries(
        [LineString([p0, p1]) for p0, p1 in zip(links["geometry0"], links["geometry1"])]
    )
    links = gpd.GeoDataFrame(links, geometry="geometry", crs=geo_crs)

    links = (
        links.assign(
            link_id=lambda df: df["bus0"] + "-" + df["bus1"] + "-DC",
            voltage=380,  # TODO Improve assumption
            length=lambda df: df.geometry.to_crs(distance_crs).length,
            underground="t",
            under_construction="f",
            tags=lambda df: df["bus0"] + " -> " + df["bus1"],
        )
        .groupby(by="link_id")
        .agg(
            {
                **{col: "first" for col in LINKS_COLUMNS if col != "p_nom"},
                "p_nom": "sum",
            }
        )[LINKS_COLUMNS]
    )
    links = gpd.GeoDataFrame(links, geometry="geometry", crs=geo_crs)

    return links


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_tyndp_network")

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    # Build node coordinates
    bidding_shapes, country_shapes = build_shapes(snakemake.input.bidding_shapes)
    buses, buses_h2 = build_buses(snakemake.input.buses, bidding_shapes, country_shapes)

    # Build links
    links = build_links(snakemake.input.reference_grid, buses)

    # Build placeholder lines, converters and transformers as empty DataFrames
    lines = gpd.GeoDataFrame(columns=LINES_COLUMNS, geometry="geometry").set_index(
        pd.Index([], name="line_id")
    )
    converters = gpd.GeoDataFrame(
        columns=CONVERTERS_COLUMNS, geometry="geometry"
    ).set_index(pd.Index([], name="converter_id"))
    transformers = gpd.GeoDataFrame(
        columns=TRANSFORMERS_COLUMNS, geometry="geometry"
    ).set_index(pd.Index([], name="transformer_id"))

    # Export to csv for base_network
    buses.to_csv(snakemake.output["substations"], quotechar="'")
    buses_h2.to_csv(snakemake.output["substations_h2"], quotechar="'")
    lines.to_csv(snakemake.output["lines"], quotechar="'")
    links.to_csv(snakemake.output["links"], quotechar="'")
    converters.to_csv(snakemake.output["converters"], quotechar="'")
    transformers.to_csv(snakemake.output["transformers"], quotechar="'")

    # Export to GeoJSON for quick validations
    buses.to_file(snakemake.output["substations_geojson"])
    buses_h2.to_file(snakemake.output["substations_h2_geojson"])
    lines.to_file(snakemake.output["lines_geojson"])
    links.to_file(snakemake.output["links_geojson"])
    converters.to_file(snakemake.output["converters_geojson"])
    transformers.to_file(snakemake.output["transformers_geojson"])
