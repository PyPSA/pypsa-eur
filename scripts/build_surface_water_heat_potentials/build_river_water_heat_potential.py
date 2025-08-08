# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Calculate river water heat potential for district heating systems.

This script computes the thermal potential of rivers as a heat source for district
heating applications. It uses HERA river discharge and ambient temperature data to
estimate available heating power and average water temperatures across regions
intersected with district heating areas.

The approximation accounts for seasonal variations in river flow and temperature,
providing both spatial and temporal aggregates. Temporal aggregates are only used for plotting.

Relevant Settings
-----------------

.. code:: yaml

    sector:
        district_heating:
            dh_area_buffer: 1000  # Buffer around DH areas in meters
            heat_source_cooling: 3  # Temperature reduction in heating system
    snapshots:
        start: "2019-01-01"
        end: "2019-12-31"
    enable:
        drop_leap_day: true

Inputs
------
- `data/hera_{year}/river_discharge_{year}.nc`: River discharge data from HERA
- `data/hera_{year}/ambient_temp_{year}.nc`: Ambient temperature data from HERA
- `resources/<run_name>/regions_onshore_base_s_{clusters}.geojson`: Onshore regions
- `resources/<run_name>/dh_areas_base_s_{clusters}.geojson`: District heating areas

Outputs
-------
- `resources/<run_name>/heat_source_power_river_water_base_s_{clusters}.csv`: River heating power potentials by region
- `resources/<run_name>/temp_river_water_base_s_{clusters}.nc`: River water temperature profiles by region
- `resources/<run_name>/temp_river_water_base_s_{clusters}_temporal_aggregate.nc`: Temporal aggregated temperature data
- `resources/<run_name>/heat_source_energy_river_water_base_s_{clusters}_temporal_aggregate.nc`: Temporal aggregated energy data
"""

import gc
import logging

import dask
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
    update_config_from_wildcards,
)
from approximators.river_water_heat_approximator import RiverWaterHeatApproximator

logger = logging.getLogger(__name__)

MEMORY_SAFETY_FACTOR = 0.7  # Use 70% of available memory for Dask arrays


def _create_empty_datasets(
    snapshots: pd.DatetimeIndex, center_lon: float, center_lat: float
) -> tuple:
    """
    Create empty spatial and temporal aggregate datasets for regions without DH areas.

    When a region has no intersection with district heating areas, we still need
    to provide valid datasets with zero values to maintain consistent data structure
    across all regions. This prevents errors in downstream processing.

    Parameters
    ----------
    snapshots : pd.DatetimeIndex
        Time snapshots for the spatial aggregate
    center_lon : float
        Longitude of region center (for fallback coordinate)
    center_lat : float
        Latitude of region center (for fallback coordinate)

    Returns
    -------
    tuple
        Tuple of (spatial_aggregate, temporal_aggregate) datasets with zero values
    """
    # Create spatial aggregate with time-series of zeros
    # This represents no available river heat power/temperature over time
    spatial_aggregate = xr.Dataset(
        data_vars={
            "total_power": xr.DataArray(
                np.zeros(len(snapshots)),  # Zero power for all timestamps
                dims=["time"],
                coords={"time": snapshots},
            ),
            "average_temperature": xr.DataArray(
                np.zeros(len(snapshots)),  # Zero temperature for all timestamps
                dims=["time"],
                coords={"time": snapshots},
            ),
        }
    )

    # Create temporal aggregate with single spatial point of zeros
    # This represents no available river heat energy/temperature spatially
    temporal_aggregate = xr.Dataset(
        data_vars={
            "total_energy": xr.DataArray(
                [[0.0]],  # Zero energy at region center
                dims=["longitude", "latitude"],
                coords={"longitude": [center_lon], "latitude": [center_lat]},
            ),
            "average_temperature": xr.DataArray(
                [[0.0]],  # Zero temperature at region center
                dims=["longitude", "latitude"],
                coords={"longitude": [center_lon], "latitude": [center_lat]},
            ),
        }
    )

    return spatial_aggregate, temporal_aggregate


def get_regional_result(
    river_discharge_fn: xr.DataArray,
    ambient_temperature_fn: xr.DataArray,
    region: gpd.GeoSeries,
    dh_areas: gpd.GeoDataFrame,
    snapshots: pd.DatetimeIndex,
    generate_temporal_aggregates: bool = False,
) -> dict:
    """
    Calculate river water heat potential for a given region and district heating areas.

    Parameters
    ----------
    river_discharge_fn : xr.DataArray or str
        Path to NetCDF file or xarray DataArray containing river discharge data.
    ambient_temperature_fn : xr.DataArray or str
        Path to NetCDF file or xarray DataArray containing ambient temperature data.
    region : geopandas.GeoSeries
        Geographical region for which to compute the heat potential.
    dh_areas : geopandas.GeoDataFrame
        District heating areas to intersect with the region.
    snapshots : pd.DatetimeIndex
        Time snapshots, used only for regions without dh_areas

    Returns
    -------
    dict
        Dictionary with keys 'spatial aggregate' and 'temporal aggregate'.
        'spatial aggregate' contains total power and average temperature.
        'temporal aggregate' contains time series for energy and temperature (for analysis/plotting).
    """
    # Store original region for fallback centroid calculation
    original_region = region.copy()

    # Intersect region with district heating areas
    intersected_geometry = gpd.overlay(
        region.to_frame(),
        dh_areas,
        how="intersection",
    ).union_all()

    region.geometry = intersected_geometry

    # Handle empty geometry case (no intersection with DH areas)
    # This occurs when a region has no district heating infrastructure
    if region.geometry.is_empty.any():
        # Get the center of the original region (before intersection)
        # We use the original region to get a meaningful coordinate for the empty datasets
        # Project to EPSG:3035 for accurate centroid calculation, then back to EPSG:4326
        region_center = (
            original_region.to_crs("EPSG:3035").centroid.to_crs("EPSG:4326").iloc[0]
        )
        center_lon = region_center.x
        center_lat = region_center.y

        # Return zero-filled datasets with proper structure
        spatial_aggregate, temporal_aggregate = _create_empty_datasets(
            snapshots, center_lon, center_lat
        )

        return {
            "spatial aggregate": spatial_aggregate,
            "temporal aggregate": temporal_aggregate,
        }

    # Process region with valid DH area intersection
    # This is the main processing path for regions with district heating

    # Get bounding box for efficient data clipping
    # We use total_bounds to get the minimal rectangular area covering all geometries
    minx, miny, maxx, maxy = region.total_bounds

    # Data processing strategy:
    # 1. Load HERA discharge and temperature data
    # 2. Clip to region bounds for efficiency
    # 3. Reproject to EPSG:3035 for accurate spatial calculations
    # 4. Feed to approximator for heat potential calculation

    # Load and preprocess river discharge data from HERA dataset
    # Use context manager to ensure proper file handle cleanup
    # Use smaller chunks to reduce memory usage
    with xr.open_dataset(
        river_discharge_fn,
        chunks={
            "time": -1,
            "lat": 14,
            "lon": 4530,
        },  # Smaller chunks to reduce memory footprint
    ) as river_ds:
        river_discharge = river_ds["dis"].chunk(
            {
                "time": -1,
                "lat": "auto",
                "lon": "auto",
            }  # Rechunk to match ambient temperature
        )

        river_discharge = (
            river_discharge
            # Standardize coordinate names to match expected format
            .rename({"lat": "latitude", "lon": "longitude"})
            # Set CRS to WGS84 for geographic operations
            .rio.write_crs("EPSG:4326")
            # Clip to region bounds to reduce memory usage and processing time
            .rio.clip_box(minx, miny, maxx, maxy)
            # Reproject to European grid (EPSG:3035) for accurate area calculations
            .rio.reproject("EPSG:3035")
        )

    # Load and preprocess ambient temperature data from HERA dataset
    # Use context manager to ensure proper file handle cleanup
    # Use smaller chunks to reduce memory usage
    logging.info("Loading ambient temperature data...")
    with xr.open_dataset(
        ambient_temperature_fn,
        chunks={
            "time": -1,
            "lat": 990,
            "lon": 1510,
        },  # Smaller chunks to reduce memory footprint
    ) as temp_ds:
        ambient_temperature = temp_ds["ta6"].chunk(
            {
                "time": -1,
                "lat": "auto",
                "lon": "auto",
            }  # Rechunk to match river discharge
        )

        ambient_temperature = (
            temp_ds["ta6"]  # Extract temperature variable (Kelvin or Celsius)
            # Apply same preprocessing pipeline as discharge data
            .rename({"lat": "latitude", "lon": "longitude"})
            .rio.write_crs("EPSG:4326")
            .rio.clip_box(minx, miny, maxx, maxy)
            .rio.reproject("EPSG:3035")
        )

    # Reproject region to match data CRS for spatial calculations
    # EPSG:3035 provides accurate area/distance calculations for Europe
    region = region.to_crs("EPSG:3035")

    # Initialize the river water heat approximator with prepared data
    river_water_heat_approximator = RiverWaterHeatApproximator(
        volume_flow=river_discharge,  # River discharge [m³/s]
        ambient_temperature=ambient_temperature,  # Air temperature [K or °C]
        region=region,  # Geographic region of interest
    )

    # Calculate spatial aggregate (time series data for the region)
    # Contains total_power [MW] and average_temperature [°C] over time
    # Keep as lazy dask array - will be computed later when needed
    spatial_aggregate = river_water_heat_approximator.get_spatial_aggregate()

    # Calculate temporal aggregate only if needed (spatial distribution data for plotting)
    if generate_temporal_aggregates:
        temporal_aggregate = (
            river_water_heat_approximator.get_temporal_aggregate()
            .rio.reproject("EPSG:4326")  # Convert back to WGS84 for output consistency
            .rename({"x": "longitude", "y": "latitude"})  # Standardize coordinate names
        )
        temporal_result = temporal_aggregate.compute()
    else:
        temporal_result = None

    # Compute results immediately to free Dask arrays and enable garbage collection
    spatial_result = spatial_aggregate.compute()

    # Explicitly delete approximator and intermediate arrays
    del river_water_heat_approximator
    del river_discharge, ambient_temperature
    if generate_temporal_aggregates:
        del temporal_aggregate
    gc.collect()

    result = {
        "spatial aggregate": spatial_result,
    }

    # Only include temporal aggregate if computed
    if temporal_result is not None:
        result["temporal aggregate"] = temporal_result

    return result


def set_dask_chunk_size(
    n_threads: int,  # Number of threads per worker,
    memory_mb: int,  # Memory per worker in MB
    memory_safety_factor=MEMORY_SAFETY_FACTOR,
    n_datasets: int = 2,  # ambient temperature and river discharge datasets
    operation_multiplier: int = 3,  # Multiplier for operation overhead
) -> None:
    """
    Set the Dask chunk size based on available memory and number of threads.
    This function calculates the chunk size for Dask arrays to optimize memory usage
    """

    chunk_size = (
        memory_mb * memory_safety_factor / n_threads / n_datasets / operation_multiplier
    )
    dask.config.set({"array.chunk-size": f"{chunk_size}MB"})


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_brownfield",
            clusters="39",
            opts="",
            ll="vopt",
            sector_opts="",
            planning_horizons=2050,
        )

    # Configure logging and scenario
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Get simulation snapshots
    snapshots: pd.DatetimeIndex = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    # Load regions and district heating areas
    regions_onshore = gpd.read_file(snakemake.input["regions_onshore"])
    regions_onshore.set_index("name", inplace=True)
    regions_onshore = regions_onshore.to_crs("EPSG:4326")

    dh_areas = gpd.read_file(snakemake.input["dh_areas"]).to_crs("EPSG:3035")
    # Buffer district heating areas by specified amount
    dh_areas["geometry"] = dh_areas.geometry.buffer(snakemake.params.dh_area_buffer)
    dh_areas = dh_areas.to_crs("EPSG:4326")

    # Configure Dask for multi-threading within operations (no distributed cluster)
    dask.config.set(scheduler="threads")  # Use threaded scheduler
    dask.config.set(num_workers=snakemake.threads)  # Use specified number of threads

    set_dask_chunk_size(
        n_threads=snakemake.threads, memory_mb=snakemake.resources.mem_mb
    )

    # Process regions sequentially but with multi-threaded Dask operations
    results = []
    for i, region_name in enumerate(regions_onshore.index, 1):
        # Extract region geometry and create a copy to avoid modification conflicts
        region = gpd.GeoSeries(regions_onshore.loc[region_name].copy(deep=True))

        # Process region with multi-threaded Dask operations
        result = get_regional_result(
            river_discharge_fn=snakemake.input.hera_river_discharge,
            ambient_temperature_fn=snakemake.input.hera_ambient_temperature,
            region=region,
            dh_areas=dh_areas,
            snapshots=snapshots,
            generate_temporal_aggregates=snakemake.params.generate_temporal_aggregates,
        )
        results.append(result)

        # Explicit cleanup to free memory between regions
        del result, region
        gc.collect()

    # Build DataFrame of total power for each region
    # This creates a time-series matrix with regions as columns and time as rows
    power = pd.DataFrame(
        {
            region_name: res["spatial aggregate"]["total_power"].to_pandas()
            for region_name, res in zip(regions_onshore.index, results)
        }
    ).dropna()  # Remove any NaN values that could cause issues downstream

    # Align power DataFrame to simulation snapshots
    # This ensures the time index exactly matches the energy system model requirements
    # Use "nearest" method to handle any minor timestamp differences due to floating point precision
    power = power.reindex(snapshots, method="nearest")

    # Save power potentials as CSV for integration with PyPSA energy system model
    # Units: MW (megawatts)
    power.to_csv(snakemake.output.heat_source_power)

    # Concatenate average temperature for all regions into single dataset
    # This creates a 2D array: [time, regions] with temperature values
    temperature = (
        xr.concat(
            [res["spatial aggregate"]["average_temperature"] for res in results],
            dim="name",
        )
        .assign_coords(name=regions_onshore.index)
        .dropna(dim="time")
    )  # Remove invalid time points

    # Align temperature data to simulation snapshots
    # This ensures temporal consistency with the energy system model
    temperature = temperature.sel(time=snapshots, method="nearest").assign_coords(
        time=snapshots
    )

    # Save temperature profiles as NetCDF for heat pump COP calculations
    # Units: °C (degrees Celsius)
    temperature.to_netcdf(snakemake.output.heat_source_temperature)

    # Save temporal aggregate results for analysis and visualization (if enabled)
    if snakemake.params.generate_temporal_aggregates:
        # Energy temporal aggregate: spatial distribution of available energy
        # Units: MWh (megawatt-hours) - total energy potential per location
        xr.concat(
            [res["temporal aggregate"]["total_energy"] for res in results],
            dim=regions_onshore.index,
        ).to_netcdf(snakemake.output.heat_source_energy_temporal_aggregate)

        # Temperature temporal aggregate: spatial distribution of temperatures
        # Units: °C (degrees Celsius) - average temperature per location
        xr.concat(
            [res["temporal aggregate"]["average_temperature"] for res in results],
            dim=regions_onshore.index,
        ).to_netcdf(snakemake.output.heat_source_temperature_temporal_aggregate)

    else:
        # Create empty placeholder files to satisfy Snakemake outputs
        # These are minimal NetCDF files that take almost no disk space
        empty_ds = xr.Dataset()
        empty_ds.to_netcdf(snakemake.output.heat_source_energy_temporal_aggregate)
        empty_ds.to_netcdf(snakemake.output.heat_source_temperature_temporal_aggregate)
