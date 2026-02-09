# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Selects n largest district heating systems from dh_areas, creates extended
onshore regions, and outputs subnode metadata for downstream rules.
"""

import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd

from scripts._helpers import configure_logging, set_scenario_config

logger = logging.getLogger(__name__)


def _sanitize_label(label):
    return (
        pd.Series([label])
        .fillna("DH")
        .astype(str)
        .str.replace(r"[^\w\s-]", "", regex=True)
        .str.strip()
        .str.replace(r"\s+", "_", regex=True)
    ).iloc[0]


def _find_containing_region(point, regions):
    """
    Find the region containing a point, or the nearest region if not contained.

    Parameters
    ----------
    point : shapely.geometry.Point
        The point to locate
    regions : geopandas.GeoDataFrame
        GeoDataFrame of regions with geometry column, indexed by region name

    Returns
    -------
    str
        Name (index) of the containing or nearest region
    """
    containing = regions[regions.contains(point)]
    if len(containing) > 0:
        return containing.index[0]
    distances = regions.geometry.distance(point)
    return distances.idxmin()


def identify_largest_district_heating_systems(
    dh_areas: gpd.GeoDataFrame,
    regions_onshore: gpd.GeoDataFrame,
    n_subnodes: int,
    countries: list[str],
    handle_missing_countries: str = "ignore",
    demand_column: str = "Dem_GWh",
    label_column: str = "Label",
) -> gpd.GeoDataFrame:
    """
    Selects the largest DH systems from dh_areas based on demand, filters by
    country, and maps each system to its containing onshore region. If a
    'city' column exists in dh_areas it is preserved for easier reporting.

    Parameters
    ----------
    dh_areas : geopandas.GeoDataFrame
        District heating areas with geometry and demand data
    regions_onshore : geopandas.GeoDataFrame
        Onshore regions indexed by cluster name
    n_subnodes : int
        Number of largest DH systems to select
    countries : list[str]
        List of country codes to filter DH areas
    handle_missing_countries : {"ignore","fill","raise"}, optional
        Behaviour if requested countries are missing in dh_areas. "ignore" logs a
        warning and proceeds with available data; "fill" currently falls back to
        ignoring (placeholder for future enrichment); "raise" stops with an error.
    demand_column : str, optional
        Column name containing demand values in GWh/a, by default "Dem_GWh"
    label_column : str, optional
        Column name containing DH system labels, by default "Label"

    Returns
    -------
    geopandas.GeoDataFrame
        Selected subnodes with columns: name, cluster, subnode_label,
        yearly_heat_demand_MWh, country, geometry, and optional city
    """
    if "country" in dh_areas.columns:
        present_countries = set(dh_areas["country"].unique())
        missing = set(countries) - present_countries
        if missing:
            msg = f"DH areas missing for countries {sorted(missing)}"
            if handle_missing_countries == "raise":
                raise ValueError(msg)
            logger.warning(msg)
        # Placeholder: if future fill logic exists, hook in here.
        dh_areas = dh_areas[
            dh_areas["country"].isin(present_countries & set(countries))
        ].copy()
        logger.info(
            f"Filtered to {len(dh_areas)} DH areas in {present_countries & set(countries)}"
        )
    else:
        logger.warning("No 'country' column in dh_areas")

    dh_areas_valid = dh_areas.dropna(subset=[demand_column])
    if len(dh_areas_valid) == 0:
        logger.warning("No valid DH areas found")
        return gpd.GeoDataFrame()

    subnodes = dh_areas_valid.nlargest(n_subnodes, demand_column).copy()
    logger.info(
        f"Selected {len(subnodes)} largest DH systems ({subnodes[demand_column].sum():.1f} GWh/a)"
    )

    subnodes["yearly_heat_demand_MWh"] = subnodes[demand_column] * 1e3

    if subnodes.crs != regions_onshore.crs:
        subnodes = subnodes.to_crs(regions_onshore.crs)

    subnodes["centroid"] = subnodes.geometry.centroid
    subnodes["cluster"] = subnodes["centroid"].apply(
        lambda p: _find_containing_region(p, regions_onshore)
    )

    subnodes["subnode_label"] = subnodes[label_column].apply(_sanitize_label)
    subnodes["name"] = subnodes["cluster"] + " " + subnodes["subnode_label"]

    keep_cols = [
        "name",
        "cluster",
        "subnode_label",
        "yearly_heat_demand_MWh",
        "country",
        "geometry",
    ]
    if "city" in subnodes.columns:
        keep_cols.insert(3, "city")

    return subnodes[keep_cols].copy()


def extend_regions_onshore(
    regions_onshore: gpd.GeoDataFrame,
    subnodes: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Extend onshore regions with subnode geometries cut from parent clusters.

    For each subnode, cuts its geometry from the parent cluster and adds
    the subnode as a new region entry.

    Parameters
    ----------
    regions_onshore : geopandas.GeoDataFrame
        Original onshore regions indexed by cluster name
    subnodes : geopandas.GeoDataFrame
        Subnodes with name, cluster, and geometry columns

    Returns
    -------
    geopandas.GeoDataFrame
        Extended regions with subnodes added and parent geometries updated
    """
    if len(subnodes) == 0:
        return regions_onshore.copy()

    subnodes_crs = subnodes.to_crs(regions_onshore.crs)
    regions_extended = regions_onshore.copy()

    for cluster in subnodes_crs["cluster"].unique():
        cluster_subnodes = subnodes_crs[subnodes_crs["cluster"] == cluster]
        subnode_union = cluster_subnodes.union_all()

        if cluster in regions_extended.index:
            original_geom = regions_extended.loc[cluster, "geometry"]
            regions_extended.loc[cluster, "geometry"] = original_geom.difference(
                subnode_union
            )

    subnode_entries = subnodes_crs[["name", "geometry"]].set_index("name")
    return pd.concat([regions_extended, subnode_entries])


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "identify_district_heating_subnodes",
            clusters=48,
        )

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    countries = snakemake.params.get("countries", [])
    subnode_countries = snakemake.params.get("subnode_countries", None)
    n_subnodes = snakemake.params.get("n_subnodes", 40)
    demand_column = snakemake.params.get("demand_column", "Dem_GWh")
    label_column = snakemake.params.get("label_column", "Label")
    handle_missing_countries = (
        snakemake.params.get(
            "handle_missing_countries",
            snakemake.config.get("sector", {})
            .get("district_heating", {})
            .get("dh_areas", {})
            .get("handle_missing_countries", "ignore"),
        )
        or "ignore"
    )

    if subnode_countries:
        invalid = set(subnode_countries) - set(countries)
        if invalid:
            logger.warning(f"Invalid subnode_countries {invalid} ignored")
        effective_subnode_countries = [c for c in subnode_countries if c in countries]
    else:
        effective_subnode_countries = countries

    logger.info(f"Identifying {n_subnodes} subnodes from {effective_subnode_countries}")

    dh_areas = gpd.read_file(snakemake.input.dh_areas)
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore).set_index("name")
    city_lookup_path = snakemake.input.get("dh_city_lookup") or snakemake.params.get(
        "dh_city_lookup", ""
    )
    if not city_lookup_path or not Path(city_lookup_path).exists():
        logger.warning(
            "dh_city_lookup missing; proceeding without city names. "
            "Run map_dh_systems_to_cities to enable city-labelled subnodes."
        )
        city_lookup = None
    else:
        city_lookup = pd.read_csv(city_lookup_path)

    subnodes = identify_largest_district_heating_systems(
        dh_areas,
        regions_onshore,
        n_subnodes,
        effective_subnode_countries,
        handle_missing_countries,
        demand_column,
        label_column,
    )

    # Attach city names from static lookup
    # Merge on both subnode_label AND country since labels are only unique within a country
    if city_lookup is not None:
        city_lookup["subnode_label"] = city_lookup["subnode_label"].apply(
            _sanitize_label
        )
        subnodes = subnodes.merge(
            city_lookup[["subnode_label", "country", "city"]],
            how="left",
            on=["subnode_label", "country"],
        )

    # Prefer city name in subnode identifier when available
    # Use the raw city name (with spaces) for readable node names
    if "city" in subnodes.columns:
        # Clean city names: remove special chars but KEEP spaces
        city_clean = (
            subnodes["city"]
            .fillna("")
            .astype(str)
            .str.replace(r"[^\w\s-]", "", regex=True)
            .str.strip()
        )
        use_city = city_clean != ""
        subnodes.loc[use_city, "name"] = (
            subnodes.loc[use_city, "cluster"] + " " + city_clean[use_city]
        )

    # Handle duplicate names
    duplicates = subnodes["name"].duplicated(keep=False)
    if duplicates.any():
        for name in subnodes.loc[duplicates, "name"].unique():
            mask = subnodes["name"] == name
            subnodes.loc[mask, "name"] = [f"{name}_{i}" for i in range(mask.sum())]

    if len(subnodes) == 0:
        logger.warning("No subnodes identified")
        subnodes = gpd.GeoDataFrame(
            columns=[
                "name",
                "cluster",
                "subnode_label",
                "city",
                "yearly_heat_demand_MWh",
                "country",
                "geometry",
            ],
            crs=regions_onshore.crs,
        )

    regions_extended = extend_regions_onshore(regions_onshore, subnodes)

    subnodes.to_file(snakemake.output.dh_subnodes, driver="GeoJSON")
    regions_extended.to_file(
        snakemake.output.regions_onshore_extended, driver="GeoJSON"
    )

    if len(subnodes) > 0:
        logger.info(
            f"Saved {len(subnodes)} subnodes ({subnodes['yearly_heat_demand_MWh'].sum() / 1e6:.2f} TWh/a)"
        )
