# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Create interactive maps for heat source temperature and energy data.

This script generates interactive Folium maps displaying heat source temperature
and energy potential data across European regions. It visualizes spatial distributions
of renewable heat sources like river water, sea water, and ambient air temperatures,
along with their energy potentials for district heating applications.

The script creates two types of maps:
- Temperature maps showing spatial temperature distributions (°C)
- Energy maps showing total energy potential (TWh) where available

Maps include regional boundaries with aggregated values and detailed point data
with interactive tooltips. Temperature data is averaged by region while energy
data is summed by region to show total potential.

Relevant Settings
-----------------

.. code:: yaml

    plotting:
        heat_source_map:
            temperature_cmap: "Reds"  # Colormap for temperature data
            energy_cmap: "Oranges"    # Colormap for energy data

Inputs
------
- `resources/<run_name>/regions_onshore_base_s_{clusters}.geojson`: Regional boundaries
- `resources/<run_name>/temp_{carrier}_base_s_{clusters}_temporal_aggregate.nc`: Temperature data
- `resources/<run_name>/heat_source_energy_{carrier}_base_s_{clusters}_temporal_aggregate.nc`: Energy data (optional)

Outputs
-------
- `results/<run_name>/plots/heat_source_map_{carrier}_temperature.html`: Interactive temperature map
- `results/<run_name>/plots/heat_source_map_{carrier}_energy.html`: Interactive energy potential map

Notes
-----
Uses Folium for interactive web-based mapping. Temperature values in °C,
energy values converted from MWh to TWh for display. Handles missing energy
data by creating empty placeholder maps.
"""

import logging

import folium
import geopandas as gpd
import xarray as xr
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
    longitude_name: str = "longitude",
    latitude_name: str = "latitude",
    onshore_region_name: str = "name",
    title: str | None = None,
    cmap: str = "viridis",
    aggregate_type: str = "mean",  # 'mean' for temperature, 'sum' for energy
) -> folium.Map:
    """
    Create an interactive Folium map from heat source data.

    Generates a web-based interactive map displaying heat source temperature or
    energy data with regional boundaries and point-based data visualization.
    Regional values are aggregated and displayed in tooltips while individual
    data points show detailed values.

    Parameters
    ----------
    da : xr.DataArray
        DataArray containing heat source data with spatial coordinates.
    regions_onshore : gpd.GeoDataFrame
        GeoDataFrame with onshore region geometries for boundary overlay.
    var_name : str
        Name of the variable to plot from the DataArray.
    longitude_name : str, default 'longitude'
        Name of the longitude coordinate in the DataArray.
    latitude_name : str, default 'latitude'
        Name of the latitude coordinate in the DataArray.
    onshore_region_name : str, default 'name'
        Column name in regions_onshore containing region identifiers.
    title : str, optional
        Title for the map legend. If None, uses var_name.
    cmap : str, default 'viridis'
        Matplotlib colormap name for data visualization.
    aggregate_type : str, default 'mean'
        Aggregation method for regional values. Use 'mean' for temperature
        data and 'sum' for energy data.

    Returns
    -------
    folium.Map
        Interactive map with layered heat source data and regional boundaries.

    Raises
    ------
    ValueError
        If required columns/coordinates are missing or invalid aggregate_type.
    """
    # Reset index if needed and check for required column
    if hasattr(regions_onshore, "index"):
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

    # Filter out zero values and NaNs
    df = df[df[var_name] != 0].dropna(subset=[var_name])
    # Convert to GeoDataFrame
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df[longitude_name], df[latitude_name]),
        crs="EPSG:4326",
    )

    # Calculate aggregations by region
    if aggregate_type == "mean":
        region_totals = gdf.groupby(onshore_region_name)[var_name].mean()
        agg_title = f"Average {var_name}"
    elif aggregate_type == "sum":  # sum for energy
        region_totals = gdf.groupby(onshore_region_name)[var_name].sum()
        agg_title = f"Total {var_name}"
    else:
        raise ValueError(
            f"Invalid aggregate_type '{aggregate_type}'. Use 'mean' or 'sum'."
        )

    # Merge region totals back to regions
    regions_with_totals = regions_onshore.merge(
        region_totals, on=onshore_region_name, how="left"
    )

    # Center the map on data
    center = [gdf[latitude_name].mean(), gdf[longitude_name].mean()]
    m = folium.Map(location=center, zoom_start=5)

    # Add region boundaries with aggregated values in the tooltip
    folium.GeoJson(
        regions_with_totals,
        name="Regions",
        style_function=lambda x: {
            "fillColor": "transparent",
            "color": "black",
            "weight": 2,
        },
        tooltip=folium.GeoJsonTooltip(
            fields=[onshore_region_name, var_name],
            aliases=[onshore_region_name, agg_title],
        ),
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
    region_onshore = regions_onshore.to_crs("EPSG:4326")

    # Get colormaps from config
    temperature_cmap = snakemake.params.plotting.get("heat_source_map", {}).get(
        "temperature_cmap", "Reds"
    )
    energy_cmap = snakemake.params.plotting.get("heat_source_map", {}).get(
        "energy_cmap", "Oranges"
    )

    logger.info(
        f"Creating temperature map for {snakemake.wildcards.carrier} heat source"
    )
    # Load temperature data
    temperature_data = xr.open_dataset(snakemake.input.heat_source_temperature)

    if snakemake.wildcards.carrier == "ambient_air":
        temp_var = "temperature"
    else:
        temp_var = "average_temperature"

    try:
        # If time dimension exists, use the mean across time
        if "time" in temperature_data[temp_var].dims:
            plot_data = temperature_data[temp_var].mean("time")
            logger.info("Taking mean across time dimension")
        else:
            plot_data = temperature_data[temp_var]

        # Create and save the temperature map
        temp_map = plot_heat_source_map(
            da=plot_data,
            regions_onshore=regions_onshore,
            var_name=temp_var,
            title=f"{snakemake.wildcards.carrier.replace('_', ' ').title()} Temperature (°C)",
            cmap=temperature_cmap,
            aggregate_type="mean",
        )

        temp_map.save(snakemake.output.temp_map)
        logger.info(f"Temperature map saved to {snakemake.output.temp_map}")

    except ValueError as e:
        logger.error(f"Error creating temperature map: {e}")

    # Plot energy map if input is available
    if (
        hasattr(snakemake.input, "heat_source_energy")
        and snakemake.input.heat_source_energy
    ):
        logger.info(
            f"Creating energy map for {snakemake.wildcards.carrier} heat source"
        )
        try:
            energy_data = xr.open_dataset(snakemake.input.heat_source_energy)

            # Create and save the energy map
            energy_map = plot_heat_source_map(
                da=energy_data["total_energy"] / 1e6,  # Convert to TWh
                regions_onshore=regions_onshore,
                var_name="total_energy",
                title=f"{snakemake.wildcards.carrier.replace('_', ' ').title()} Energy Potential (TWh)",
                cmap=energy_cmap,
                aggregate_type="sum",
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
