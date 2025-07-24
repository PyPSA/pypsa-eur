# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Calculate sea water heat potential for district heating systems.

This script computes the thermal potential of sea water as a heat source for district
heating applications. It uses sea water temperature data to estimate average water
temperatures across regions intersected with district heating areas.

The approximation provides temperature profiles for coastal regions with access to
sea water, supporting both spatial and temporal aggregates. Temporal aggregates are only used for plotting.

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
- `data/seawater_temperature.nc`: Sea water temperature data
- `resources/<run_name>/regions_onshore_base_s_{clusters}.geojson`: Onshore regions
- `resources/<run_name>/dh_areas_base_s_{clusters}.geojson`: District heating areas

Outputs
-------
- `resources/<run_name>/temp_sea_water_base_s_{clusters}.nc`: Sea water temperature profiles by region
- `resources/<run_name>/temp_sea_water_base_s_{clusters}_temporal_aggregate.nc`: Temporal aggregated temperature data
"""

import logging

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
from approximators.sea_water_heat_approximator import SeaWaterHeatApproximator
from dask.distributed import Client, LocalCluster

logger = logging.getLogger(__name__)


def _create_empty_datasets(
    snapshots: pd.DatetimeIndex, center_lon: float, center_lat: float
) -> tuple:
    """
    Create empty spatial and temporal aggregate datasets for regions without DH areas.

    When a region has no intersection with district heating areas (e.g., landlocked
    regions for sea water), we still need to provide valid datasets with zero values
    to maintain consistent data structure across all regions.

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
    # Represents no available sea water temperature over time
    spatial_aggregate = xr.Dataset(
        data_vars={
            "average_temperature": xr.DataArray(
                np.zeros(len(snapshots)),  # Zero temperature for all timestamps
                dims=["time"],
                coords={"time": snapshots},
            ),
        }
    )

    # Create temporal aggregate with single spatial point of zeros
    # Represents no available sea water temperature spatially
    temporal_aggregate = xr.Dataset(
        data_vars={
            "average_temperature": xr.DataArray(
                [[0.0]],  # Zero temperature at region center
                dims=["longitude", "latitude"],
                coords={"longitude": [center_lon], "latitude": [center_lat]},
            ),
        }
    )

    return spatial_aggregate, temporal_aggregate


def get_regional_result(
    seawater_temperature_fn: str,
    region: gpd.GeoSeries,
    dh_areas: gpd.GeoDataFrame,
    snapshots: pd.DatetimeIndex,
) -> dict:
    """
    Calculate sea water heat potential for a given region and district heating areas.

    Parameters
    ----------
    seawater_temperature_fn : str
        Path to NetCDF file containing sea water temperature data.
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
        'spatial aggregate' contains average temperature.
        'temporal aggregate' contains time series for temperature (for analysis/plotting).
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
    # This occurs when coastal regions have no district heating or inland regions
    # are processed (which shouldn't have sea water access)
    if region.geometry.is_empty.any():
        logger.info(
            "Region has no DH areas or no sea access - returning empty datasets"
        )

        # Get the center of the original region (before intersection)
        # We use the original region to get a meaningful coordinate for the empty datasets
        region_center = original_region.to_crs("EPSG:4326").centroid.iloc[0]
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
    # This is the main processing path for coastal regions with district heating

    # Get bounding box for efficient data clipping
    # We use total_bounds to get the minimal rectangular area covering all geometries
    minx, miny, maxx, maxy = region.total_bounds

    # Load and preprocess sea water temperature data
    # Data processing strategy:
    # 1. Load sea water temperature from NetCDF file
    # 2. Average over depth dimension (surface temperature)
    # 3. Clip to region bounds for efficiency
    # 4. Reproject to EPSG:3035 for accurate spatial calculations
    water_temperature = (
        xr.open_dataset(
            seawater_temperature_fn,
            chunks={
                "time": "auto",  # Chunk time dimension for memory efficiency
                "latitude": "auto",  # Chunk spatial dimensions
                "longitude": "auto",
                "depth": 1,  # Single depth chunk (usually surface only)
            },
            # engine="rasterio",  # Alternative engine for some data formats
        )["thetao"]  # Extract sea water potential temperature
        .mean(dim="depth")  # Average over depth to get surface temperature
        .rio.write_crs("EPSG:4326")  # Set CRS to WGS84 for geographic operations
        .rio.clip_box(minx, miny, maxx, maxy)  # Clip to region bounds
        .rio.reproject("EPSG:3035")  # Reproject to European grid for area calculations
    )

    # Reproject region to match data CRS for spatial calculations
    # EPSG:3035 provides accurate area/distance calculations for Europe
    region = region.to_crs("EPSG:3035")

    # Initialize the sea water heat approximator with prepared data
    # This class handles the temperature analysis for coastal heat pumps
    seawater_heat_approximator = SeaWaterHeatApproximator(
        water_temperature=water_temperature,  # Sea water temperature [°C or K]
        region=region,  # Geographic region of interest
    )

    return {
        # Calculate spatial aggregate (time series data for the region)
        # Contains average_temperature [°C] over time for heat pump COP calculations
        "spatial aggregate": seawater_heat_approximator.get_spatial_aggregate().compute(),
        # Calculate temporal aggregate (spatial distribution data)
        # Contains temperature maps for analysis/plotting - converted back to WGS84
        "temporal aggregate": seawater_heat_approximator.get_temporal_aggregate()
        .rio.reproject("EPSG:4326")  # Convert back to WGS84 for output consistency
        .rename({"x": "longitude", "y": "latitude"})  # Standardize coordinate names
        .compute(),  # Execute all lazy operations
    }


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

    # Configure logging and scenario settings from Snakemake
    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    # Get simulation time snapshots from configuration
    # This defines the temporal resolution of the energy system model
    snapshots: pd.DatetimeIndex = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    # Load geographic data for processing
    # Load onshore regions (countries/NUTS regions) for sea water analysis
    regions_onshore = gpd.read_file(snakemake.input["regions_onshore"])
    regions_onshore.set_index("name", inplace=True)  # Use region name as index
    regions_onshore.set_crs("EPSG:4326", inplace=True)  # Ensure WGS84 CRS

    # Load and preprocess district heating areas
    # These define where sea water heat pumps could be connected
    dh_areas = gpd.read_file(snakemake.input["dh_areas"]).to_crs("EPSG:3035")
    # Buffer DH areas to account for reasonable transport distances
    # This allows sea water access even if DH areas don't directly touch the coast
    dh_areas["geometry"] = dh_areas.geometry.buffer(snakemake.params.dh_area_buffer)
    dh_areas = dh_areas.to_crs("EPSG:4326")  # Convert back to WGS84 for intersection

    # Set up Dask cluster for parallel computation across regions
    # This enables processing multiple regions simultaneously for better performance
    cluster = LocalCluster(
        n_workers=int(snakemake.threads),  # One worker per available thread
        threads_per_worker=1,  # Single-threaded workers to avoid conflicts
        memory_limit=f"{snakemake.resources.mem_mb / snakemake.threads}MB",  # Distribute memory evenly
    )
    client = Client(cluster)

    # Process each region in parallel using Dask distributed computing
    # Each region is processed independently to calculate its sea water temperature
    futures = []
    for region_name in regions_onshore.index:
        logging.info(f"Processing region {region_name}")

        # Extract region geometry and create a copy to avoid modification conflicts
        region = gpd.GeoSeries(regions_onshore.loc[region_name].copy(deep=True))

        # Submit region processing task to Dask cluster
        # Each task will:
        # 1. Intersect region with DH areas (coastal access check)
        # 2. Load and clip sea water temperature data to region bounds
        # 3. Calculate temperature profiles using approximator
        futures.append(
            get_regional_result(
                seawater_temperature_fn=snakemake.input["seawater_temperature"],
                region=region,
                dh_areas=dh_areas,
                snapshots=snapshots,
            )
        )

    # Collect results from all parallel computations
    # This blocks until all regions are processed
    results = client.gather(futures)

    # Concatenate average temperature for all regions into single dataset
    # This creates a 2D array: [time, regions] with sea water temperature values
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
    # Use "nearest" method to handle any minor timestamp differences
    temperature = temperature.sel(time=snapshots, method="nearest").assign_coords(
        time=snapshots
    )

    # Save temperature profiles as NetCDF for heat pump COP calculations
    # Units: °C (degrees Celsius) - sea water temperature for coastal heat pumps
    temperature.to_netcdf(snakemake.output.heat_source_temperature)

    # Save temporal aggregate results for analysis and visualization
    # This is a spatial map showing temperature distribution along coastlines

    # Temperature temporal aggregate: spatial distribution of sea water temperatures
    # Units: °C (degrees Celsius) - average temperature per coastal location
    # Used for detailed analysis and plotting of coastal heat pump potential
    xr.concat(
        [res["temporal aggregate"]["average_temperature"] for res in results],
        dim=regions_onshore.index,
    ).to_netcdf(snakemake.output.heat_source_temperature_temporal_aggregate)
