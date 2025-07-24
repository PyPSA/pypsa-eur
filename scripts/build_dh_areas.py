# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""Build district heating areas with missing country handling.

This script processes district heating areas data and handles missing countries
according to the configured strategy (ignore, fill, or raise).
"""

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon
from scripts._helpers import set_scenario_config



def handle_missing_countries(dh_areas, regions_onshore, missing_countries, handle_mode):
    """Add or fill missing countries in dh_areas according to handle_mode.
    
    Parameters
    ----------
    dh_areas : gpd.GeoDataFrame
        Existing district heating areas data.
    regions_onshore : gpd.GeoDataFrame
        Onshore regions data to use for filling missing countries.
    missing_countries : pd.Index
        Index of country codes that are missing from dh_areas.
    handle_mode : str
        How to handle missing countries. Options:
        - 'ignore': Add countries with no geometry
        - 'fill': Add countries with onshore region geometry  
        - 'raise': Raise an error
        
    Returns
    -------
    gpd.GeoDataFrame
        Updated dh_areas with missing countries handled.
        
    Raises
    ------
    ValueError
        If handle_mode is 'raise' or an invalid mode is provided.
    """
    # Create new rows based on handling mode
    if handle_mode == "ignore":
        return dh_areas # No changes, just return original GeoDataFrame
    elif handle_mode == "fill":
        # Add countries using their onshore region geometry
        new_rows = []
        for country_code in missing_countries:
            # Filter regions for this country
            country_regions = regions_onshore[regions_onshore['name'].str.startswith(country_code)]
            country_geometry = country_regions.union_all()
            
            # Convert geometry to match dh_areas CRS
            if country_geometry is not None and not country_geometry.is_empty:
                # Create a temporary GeoSeries to handle CRS conversion
                temp_geoseries = gpd.GeoSeries([country_geometry], crs=regions_onshore.crs)
                country_geometry = temp_geoseries.to_crs(dh_areas.crs).iloc[0]
            
            new_rows.append({
                "Label": np.nan,
                "country": country_code,
                "Dem_GWh": np.nan,
                "geometry": country_geometry
            })

        # Add new rows to the GeoDataFrame if any were created
        new_gdf = gpd.GeoDataFrame(new_rows, geometry="geometry", crs=dh_areas.crs)
        return pd.concat([dh_areas, new_gdf], ignore_index=True)
    elif handle_mode == "raise":
        raise ValueError(
            f"The following countries are missing in the district heating areas data: {missing_countries.to_list()}. Set `config:sector:district_heating:   handle_missing_countries` to 'ignore' to assume no district heating areas in these countries or 'fill' if you want to assume they have the same district heating areas as their onshore regions. `dh_areas` are used for the computation of some heat sources and storage potentials but not heat load."
        )
    else:
        raise ValueError(
            f"Invalid value for `config:sector:district_heating:   handle_missing_countries`: {handle_mode}. Valid values are 'ignore', 'fill', or 'raise'."
        )

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            clusters=48,
            planning_horizons="2050",
        )

    set_scenario_config(snakemake)

    # Load input data
    dh_areas: gpd.GeoDataFrame = gpd.read_file(snakemake.input.dh_areas)
    regions_onshore: gpd.GeoDataFrame = gpd.read_file(snakemake.input.regions_onshore)

    # Identify countries present in onshore regions but missing from dh_areas
    region_countries = set([name.split()[0][:2] for name in regions_onshore['name']])
    dh_countries = set(dh_areas['country'].unique())
    missing_countries = pd.Index(list(region_countries - dh_countries))


    # Handle missing countries according to configuration
    if not missing_countries.empty:
        dh_areas = handle_missing_countries(
            dh_areas,
            regions_onshore,
            missing_countries,
            snakemake.params["handle_missing_countries"]
        )

    # Save the processed district heating areas
    dh_areas.to_file(snakemake.output.dh_areas, driver="GeoJSON")
