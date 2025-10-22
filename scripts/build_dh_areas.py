# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Build and validate district heating areas for energy system modeling.
District heating areas are used for computing heat source and storage potentials.

This script processes district heating (DH) areas data from external sources and
ensures all modeled countries have consistent representation. It handles missing
countries that exist in the onshore regions but lack district heating data according
to configurable strategies.

The script supports three strategies for handling missing countries:
- 'ignore': Countries without DH data are assumed to have no district heating
- 'fill': Countries are assigned their full onshore region as potential DH area
- 'raise': Missing countries cause an error to ensure explicit handling

Relevant Settings
-----------------

.. code:: yaml

    countries: ['DE', 'FR', 'ES', ...]  # List of modeled countries
    sector:
        district_heating:
            dh_areas:
                handle_missing_countries: 'ignore'  # or 'fill' or 'raise'

Inputs
------
- `data/dh_areas.gpkg`: District heating areas data (GeoPackage format)
- `resources/<run_name>/regions_onshore_base_s_{clusters}.geojson`: Onshore regions for reference

Outputs
-------
- `resources/<run_name>/dh_areas_base_s_{clusters}.geojson`: Processed district heating areas with missing countries handled
"""

import logging

import geopandas as gpd
import numpy as np
import pandas as pd

from scripts._helpers import set_scenario_config

logger = logging.getLogger(__name__)


def handle_missing_countries(dh_areas, regions_onshore, missing_countries, handle_mode):
    """
    Handle countries that exist in onshore regions but lack district heating data.

    This function ensures consistency between the modeled countries and the available
    district heating areas data. It's common for some countries to lack detailed DH
    infrastructure data, which needs to be handled explicitly to avoid modeling errors.

    Parameters
    ----------
    dh_areas : gpd.GeoDataFrame
        Existing district heating areas data with columns: Label, country, Dem_GWh, geometry
    regions_onshore : gpd.GeoDataFrame
        Onshore regions data to use for filling missing countries
    missing_countries : pd.Index
        Index of country codes (2-letter ISO codes) missing from dh_areas
    handle_mode : str
        Strategy for handling missing countries:
        - 'ignore': Assume no district heating exists (returns unchanged data)
        - 'fill': Use full onshore region as potential DH area
        - 'raise': Fail with informative error message

    Returns
    -------
    gpd.GeoDataFrame
        Updated dh_areas with missing countries handled according to strategy

    Raises
    ------
    ValueError
        If handle_mode is 'raise' or an invalid mode is provided
    """
    if handle_mode == "ignore":
        # Strategy: Assume missing countries have no district heating infrastructure
        # This is conservative but may underestimate DH potential in some regions
        return dh_areas  # No changes needed - missing countries simply won't appear

    elif handle_mode == "fill":
        # Strategy: Use full onshore region geometry as potential DH area
        # This is optimistic but ensures no country is excluded
        # Create entries for missing countries using their onshore region boundaries
        new_rows = []
        for country_code in missing_countries:
            # Find all regions belonging to this country
            # Region names typically start with 2-letter country code (e.g., 'DE 1', 'FR 2')
            country_regions = regions_onshore[
                regions_onshore["name"].str.startswith(country_code)
            ]

            # Merge all regions of this country into a single geometry
            # This creates the maximum possible DH area for the country
            country_geometry = country_regions.union_all()

            # Handle CRS conversion using temporary GeoSeries
            if country_geometry is not None and not country_geometry.is_empty:
                temp_geoseries = gpd.GeoSeries(
                    [country_geometry], crs=regions_onshore.crs
                )
                country_geometry = temp_geoseries.to_crs(dh_areas.crs).iloc[0]

            # Create new row with country's full onshore area as potential DH area
            new_rows.append(
                {
                    "Label": np.nan,  # No specific DH area label
                    "country": country_code,  # 2-letter ISO country code
                    "Dem_GWh": np.nan,  # No demand data available
                    "geometry": country_geometry,  # Full country geometry
                }
            )

        # Append new country entries to existing DH areas data
        if new_rows:  # Only create GeoDataFrame if there are rows to add
            new_gdf = gpd.GeoDataFrame(new_rows, geometry="geometry", crs=dh_areas.crs)
            return pd.concat([dh_areas, new_gdf], ignore_index=True)
        else:
            return dh_areas  # No missing countries to add
    elif handle_mode == "raise":
        # Strategy: Fail explicitly to force manual handling of missing countries
        # This ensures users are aware of data gaps and make conscious decisions
        raise ValueError(
            f"Missing district heating data for countries: {missing_countries.to_list()}. "
            f"Configure 'sector.district_heating.dh_areas.handle_missing_countries' to:\n"
            f"  - 'ignore': Assume no DH infrastructure in missing countries\n"
            f"  - 'fill': Use full onshore regions as potential DH areas\n"
            f"  - 'raise': Current setting - requires explicit data or configuration\n\n"
            f"Note: DH areas affect heat source/storage potential calculations but not heat demand."
        )
    else:
        # Invalid configuration
        raise ValueError(
            f"Invalid handle_missing_countries setting: '{handle_mode}'. "
            f"Valid options: 'ignore', 'fill', 'raise'"
        )


if __name__ == "__main__":
    # Mock snakemake object for development/testing
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_dh_areas",  # Rule name for testing
            clusters=48,  # Number of network clusters
        )

    # Apply scenario configuration from Snakemake workflow
    set_scenario_config(snakemake)

    # Load input geographic data
    # District heating areas: contains existing DH infrastructure boundaries
    dh_areas: gpd.GeoDataFrame = gpd.read_file(snakemake.input.dh_areas)

    # Onshore regions: contains all modeled country/region boundaries
    # Used as reference for identifying missing countries and potential fill geometries
    regions_onshore: gpd.GeoDataFrame = gpd.read_file(snakemake.input.regions_onshore)

    # Identify discrepancies between modeled countries and available DH data
    # Extract country codes from region names (assumes format like 'DE 1', 'FR 2', etc.)
    region_countries = set([name.split()[0][:2] for name in regions_onshore["name"]])

    # Get countries that already have DH area data
    dh_countries = set(dh_areas["country"].unique())

    # Find countries in the model but missing from DH areas data
    # These countries need to be handled according to the configured strategy
    missing_countries = pd.Index(list(region_countries - dh_countries))

    # Process missing countries according to configured strategy
    if not missing_countries.empty:
        logger.info(
            f"Found {len(missing_countries)} missing countries: {list(missing_countries)}"
        )
        logger.info(
            f"Handling strategy: {snakemake.params['handle_missing_countries']}"
        )

        dh_areas = handle_missing_countries(
            dh_areas,
            regions_onshore,
            missing_countries,
            snakemake.params["handle_missing_countries"],
        )
    else:
        logger.info("All modeled countries have district heating areas data")

    # Save the processed district heating areas for downstream use
    # Output format: GeoJSON for compatibility with other PyPSA-Eur scripts
    dh_areas.to_file(snakemake.output.dh_areas, driver="GeoJSON")
    logger.info(
        f"Saved {len(dh_areas)} district heating areas to {snakemake.output.dh_areas}"
    )
