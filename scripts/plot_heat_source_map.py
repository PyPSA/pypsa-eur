# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Plot heat source temperature and energy maps as interactive HTML files.
"""

import logging
import os
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr
import folium
from shapely.geometry import Point
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)

logger = logging.getLogger(__name__)


def plot_heat_source_map(
    da: xr.DataArray, 
    regions_onshore: gpd.GeoDataFrame, 
    var_name: str, 
    longitude_name: str="longitude", 
    latitude_name: str="latitude", 
    onshore_region_name: str="name",
    title: str=None,
    cmap: str="viridis"
):
    """
    Create an interactive folium map from a DataArray with heat source data.
    
    Parameters
    ----------
    da : xr.DataArray
        Data array containing the heat source data
    regions_onshore : gpd.GeoDataFrame
        GeoDataFrame with onshore region geometries
    var_name : str
        Name of the variable to plot
    longitude_name : str, default 'longitude'
        Name of the longitude coordinate in the DataArray
    latitude_name : str, default 'latitude'
        Name of the latitude coordinate in the DataArray
    onshore_region_name : str, default 'name'
        Name of the column in regions_onshore containing region identifiers
    title : str, optional
        Title for the map
    cmap : str, default 'viridis'
        Colormap to use for the plot
        
    Returns
    -------
    folium.Map
        Interactive map with heat source data
    """
    # Reset index if needed and check for required column
    if hasattr(regions_onshore, 'index') and onshore_region_name == 'name':
        regions_onshore = regions_onshore.reset_index()
        
    if onshore_region_name not in regions_onshore.columns:
        raise ValueError(f"Column '{onshore_region_name}' not found in regions_onshore")
        
    # Convert DataArray to DataFrame
    df = da.to_dataframe().reset_index()
    
    # Check that required coordinate names exist
    if longitude_name not in df.columns:
        raise ValueError(f"Coordinate '{longitude_name}' not found in DataArray")
    
    if latitude_name not in df.columns:
        raise ValueError(f"Coordinate '{latitude_name}' not found in DataArray")
    
    # Check that variable name exists
    if var_name not in df.columns:
        raise ValueError(f"Variable '{var_name}' not found in DataArray")
    
    # Create a GeoDataFrame by constructing a Point for each (longitude, latitude) pair
    gdf = gpd.GeoDataFrame(
        df,
        geometry=[Point(lon, lat) for lon, lat in zip(df[longitude_name], df[latitude_name])],
        crs="EPSG:4326",
    ).replace(0, np.nan).dropna(subset=[var_name])
    
    # Center the map on data
    center = [gdf[latitude_name].mean(), gdf[longitude_name].mean()]
    m = folium.Map(location=center, zoom_start=5)
    
    # Add region boundaries
    folium.GeoJson(
        regions_onshore,
        name="Regions",
        style_function=lambda x: {"fillColor": "transparent", "color": "black", "weight": 2},
        tooltip=folium.GeoJsonTooltip(fields=[onshore_region_name]),
    ).add_to(m)
    
    # Add data points
    gdf.explore(
        m=m,
        column=var_name,
        cmap=cmap,
        legend=True,
        legend_kwds={"caption": var_name if title is None else title},
        tooltip=[var_name],
        name=var_name,
        markersize=50,
    )
    
    # Add layer control
    folium.LayerControl().add_to(m)
    
    return m


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_heat_source_map",
            clusters="39",
            opts="",
            sector_opts="",
            planning_horizons="2050",
            carrier="river_water",
        )
    
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)
    
    # Load onshore regions shapefile
    regions_onshore = gpd.read_file(snakemake.input.regions)
    
    # Get colormaps from config
    temperature_cmap = snakemake.params.plotting.get("heat_source_map", {}).get("temperature_cmap", "Reds")
    energy_cmap = snakemake.params.plotting.get("heat_source_map", {}).get("energy_cmap", "Oranges")
    
    # Plot temperature map
    logger.info(f"Creating temperature map for {snakemake.wildcards.carrier} heat source")
    temperature_data = xr.open_dataset(snakemake.input.heat_source_temperature)
    
    # Find the first data variable in the dataset
    temp_vars = list(temperature_data.data_vars)
    if not temp_vars:
        logger.error("No temperature variables found in the dataset")
    else:
        temp_var = temp_vars[0]
        logger.info(f"Found temperature variable: {temp_var}")
        
        try:
            # If time dimension exists, use the mean across time
            if 'time' in temperature_data[temp_var].dims:
                plot_data = temperature_data[temp_var].mean('time')
                logger.info("Taking mean across time dimension")
            else:
                plot_data = temperature_data[temp_var]
                
            # Create and save the temperature map
            temp_map = plot_heat_source_map(
                plot_data,
                regions_onshore,
                temp_var,
                title=f"{snakemake.wildcards.carrier.replace('_', ' ').title()} Temperature (°C)",
                cmap=temperature_cmap
            )
            
            temp_map.save(snakemake.output.temp_map)
            logger.info(f"Temperature map saved to {snakemake.output.temp_map}")
            
        except ValueError as e:
            logger.error(f"Error creating temperature map: {e}")
    
    # Plot energy map if input is available
    if hasattr(snakemake.input, 'heat_source_energy') and snakemake.input.heat_source_energy:
        logger.info(f"Creating energy map for {snakemake.wildcards.carrier} heat source")
        try:
            energy_data = xr.open_dataset(snakemake.input.heat_source_energy)
            
            # Find the first data variable in the dataset
            energy_vars = list(energy_data.data_vars)
            if not energy_vars:
                logger.error("No energy variables found in the dataset")
                # Create an empty map to satisfy output requirements
                empty_map = folium.Map(location=[50, 10], zoom_start=4)
                empty_map.save(snakemake.output.energy_map)
                logger.warning(f"Created empty energy map at {snakemake.output.energy_map}")
            else:
                energy_var = energy_vars[0]
                logger.info(f"Found energy variable: {energy_var}")
                
                # Create and save the energy map
                energy_map = plot_heat_source_map(
                    energy_data[energy_var],
                    regions_onshore,
                    energy_var,
                    title=f"{snakemake.wildcards.carrier.replace('_', ' ').title()} Energy Potential (MWh)",
                    cmap=energy_cmap
                )
                
                energy_map.save(snakemake.output.energy_map)
                logger.info(f"Energy map saved to {snakemake.output.energy_map}")
        except Exception as e:
            logger.error(f"Error creating energy map: {e}")
            # Create an empty map to satisfy output requirements
            empty_map = folium.Map(location=[50, 10], zoom_start=4)
            empty_map.save(snakemake.output.energy_map)
            logger.warning(f"Created empty energy map at {snakemake.output.energy_map}")
    else:
        # If no energy data is available, create an empty map
        logger.warning(f"No energy data available for {snakemake.wildcards.carrier}")
        empty_map = folium.Map(location=[50, 10], zoom_start=4)
        empty_map.save(snakemake.output.energy_map)
        logger.warning(f"Created empty energy map at {snakemake.output.energy_map}")
