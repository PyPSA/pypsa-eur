# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Calculate river water heat potential for district heating systems.

This script computes the thermal potential of rivers as a heat source for district
heating applications. It uses HERA river discharge and ambient temperature data to
estimate available heating power and average water temperatures across regions
intersected with district heating areas.

The approximation accounts for temporal and spatial variations in river flow and temperature,
providing both spatial and temporal aggregates. Temporal aggregates are only used for plotting.

Relevant Settings
-----------------

.. code:: yaml

    sector:
        district_heating:
            dh_area_buffer: # Buffer around DH areas in meters to include nearby rivers
            heat_source_cooling: # Exploitable temperature delta
    snapshots:
        start:
        end:
    enable:
        drop_leap_day:

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


def load_hera_data(
    hera_inputs: dict,
    snapshots: pd.DatetimeIndex,
    minx: float,
    miny: float,
    maxx: float,
    maxy: float,
) -> dict:
    """
    Load and concatenate HERA data files for multiple years with spatial clipping.

    Parameters
    ----------
    hera_inputs : dict
        Dictionary with year-specific HERA file paths from input_hera_data().
        Expected keys: hera_river_discharge_{year}, hera_ambient_temperature_{year}
    snapshots : pd.DatetimeIndex
        Target snapshots to select from the combined data
    minx, miny, maxx, maxy : float
        Bounding box coordinates for spatial clipping

    Returns
    -------
    dict
        Dictionary with processed xarray DataArrays for river_discharge and ambient_temperature
    """
    # Group files by data type
    river_files = [
        v for k, v in hera_inputs.items() if k.startswith("hera_river_discharge_")
    ]
    temp_files = [
        v for k, v in hera_inputs.items() if k.startswith("hera_ambient_temperature_")
    ]

    # Determine time range from snapshots with buffer to ensure HERA coverage
    # HERA data is on 6h intervals, so we need to extend the range to capture
    # HERA timestamps that bracket our actual snapshots
    buffer = pd.Timedelta(hours=12)  # Buffer to ensure we capture surrounding HERA data
    start_time = snapshots.min() - buffer
    end_time = snapshots.max() + buffer

    result = {}

    # Load and concatenate river discharge files using open_mfdataset
    river_discharge = xr.open_mfdataset(
        river_files,
        chunks={"time": -1, "lat": 14, "lon": 4530},
        concat_dim="time",
        combine="nested",
    )["dis"]

    # Select time range that covers our snapshots (using native HERA resolution)
    river_discharge = river_discharge.sel(time=slice(start_time, end_time))

    # Process river discharge data
    river_discharge = (
        river_discharge.rename({"lat": "latitude", "lon": "longitude"})
        .rio.write_crs("EPSG:4326")
        .rio.clip_box(minx, miny, maxx, maxy)
        .rio.reproject("EPSG:3035")
    )
    result["river_discharge"] = river_discharge

    # Load and concatenate ambient temperature files using open_mfdataset
    ambient_temperature = xr.open_mfdataset(
        temp_files,
        chunks={"time": -1, "lat": 990, "lon": 1510},
        concat_dim="time",
        combine="nested",
    )["ta6"]

    # Select time range that covers our snapshots (using native HERA resolution)
    ambient_temperature = ambient_temperature.sel(time=slice(start_time, end_time))

    # Process ambient temperature data
    ambient_temperature = (
        ambient_temperature.rename({"lat": "latitude", "lon": "longitude"})
        .rio.write_crs("EPSG:4326")
        .rio.clip_box(minx, miny, maxx, maxy)
        .rio.reproject("EPSG:3035")
    )
    result["ambient_temperature"] = ambient_temperature

    return result


def _create_empty_datasets(
    snapshots: pd.DatetimeIndex, center_lon: float, center_lat: float
) -> tuple[xr.Dataset, xr.Dataset]:
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
    tuple[xr.Dataset, xr.Dataset]
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
    hera_inputs: dict,
    region: gpd.GeoSeries,
    dh_areas: gpd.GeoDataFrame,
    snapshots: pd.DatetimeIndex,
    enable_heat_source_maps: bool = False,
) -> dict[str, xr.Dataset]:
    """
    Calculate river water heat potential for a given region and district heating areas.

    Parameters
    ----------
    hera_inputs : dict
        Dictionary containing HERA input file paths with year-specific keys.
    region : geopandas.GeoSeries
        Geographical region for which to compute the heat potential.
    dh_areas : geopandas.GeoDataFrame
        District heating areas to intersect with the region.
    snapshots : pd.DatetimeIndex
        Time snapshots, used for loading data and for regions without dh_areas
    enable_heat_source_maps : bool, default False
        Whether to enable heat source mapping.

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
    # This occurs when a region has no district heating area
    # Note: the region could still have district heating, as central heat demand is computed elsewhere
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
    # Get bounding box for efficient data clipping
    minx, miny, maxx, maxy = region.total_bounds

    # Data processing strategy:
    # 1. Load HERA discharge and temperature data
    # 2. Clip to region bounds for efficiency
    # 3. Reproject to EPSG:3035 for accurate spatial calculations
    # 4. Feed to approximator for heat potential calculation

    # Load and concatenate HERA data for all required years with preprocessing
    hera_data = load_hera_data(hera_inputs, snapshots, minx, miny, maxx, maxy)
    river_discharge = hera_data["river_discharge"]
    ambient_temperature = hera_data["ambient_temperature"]

    # Reproject region to match data CRS for spatial calculations
    region = region.to_crs("EPSG:3035")

    river_water_heat_approximator = RiverWaterHeatApproximator(
        volume_flow=river_discharge,  # River discharge (volume flow)
        ambient_temperature=ambient_temperature,  # Air temperature
        region=region,  # Geographic region of interest
    )

    # Calculate spatial aggregate (time series data for the entire region)
    # Contains total_power [MW] and average_temperature [°C] over time
    spatial_aggregate = river_water_heat_approximator.get_spatial_aggregate()

    # Calculate temporal aggregate only if needed (spatial distribution data for plotting, no time dimension)
    if enable_heat_source_maps:
        temporal_aggregate = (
            river_water_heat_approximator.get_temporal_aggregate()
            .rio.reproject("EPSG:4326")  # Convert back to WGS84 for output consistency
            .rename({"x": "longitude", "y": "latitude"})  # Standardize coordinate names
        )
        temporal_aggregate = temporal_aggregate.compute()
    else:
        temporal_aggregate = None

    # Compute results immediately to free Dask arrays and enable garbage collection
    spatial_aggregate = spatial_aggregate.compute()

    # Explicitly delete approximator and intermediate arrays
    del river_water_heat_approximator
    del river_discharge, ambient_temperature
    gc.collect()

    result = {
        "spatial aggregate": spatial_aggregate,
    }

    # Only include temporal aggregate if computed
    if temporal_aggregate is not None:
        result["temporal aggregate"] = temporal_aggregate
    if enable_heat_source_maps:
        del temporal_aggregate

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
            "build_river_water_heat_potential",
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
            hera_inputs=dict(snakemake.input),
            region=region,
            dh_areas=dh_areas,
            snapshots=snapshots,
            enable_heat_source_maps=snakemake.params.enable_heat_source_maps,
        )
        results.append(result)

        # Explicit cleanup to free memory between regions
        del result, region
        gc.collect()

    # Build DataFrame of total power for each region
    # Regions as columns and time as rows
    power = pd.DataFrame(
        {
            region_name: res["spatial aggregate"]["total_power"].to_pandas()
            for region_name, res in zip(regions_onshore.index, results)
        }
    ).dropna()

    power = power.reindex(
        snapshots, method="nearest"
    )  # Use "nearest" method to handle any minor timestamp differences due to floating point precision

    # Save power potentials in MW
    power.to_csv(snakemake.output.heat_source_power)

    # Concatenate average temperature for all regions into single dataset
    temperature = (
        xr.concat(
            [res["spatial aggregate"]["average_temperature"] for res in results],
            dim="name",
        )
        .assign_coords(name=regions_onshore.index)
        .dropna(dim="time")
    )

    # Align temperature data to snapshots, use nearest to handle any minor decimal differences
    temperature = temperature.sel(time=snapshots, method="nearest").assign_coords(
        time=snapshots
    )

    # Save temperature profiles as NetCDF for heat pump COP calculations
    # Units: °C (degrees Celsius)
    temperature.to_netcdf(snakemake.output.heat_source_temperature)

    # Save temporal aggregate results for analysis and visualization (if enabled)
    if snakemake.params.enable_heat_source_maps:
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
        empty_ds = xr.Dataset()
        empty_ds.to_netcdf(snakemake.output.heat_source_energy_temporal_aggregate)
        empty_ds.to_netcdf(snakemake.output.heat_source_temperature_temporal_aggregate)
