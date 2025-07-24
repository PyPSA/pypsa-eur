# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import logging

import geopandas as gpd
import pandas as pd
import xarray as xr
import numpy as np
from dask.distributed import Client, LocalCluster
import shapely.geometry
from _helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
    update_config_from_wildcards,
)
from approximators.river_water_heat_approximator import RiverWaterHeatApproximator

logger = logging.getLogger(__name__)


def _create_empty_datasets(snapshots: pd.DatetimeIndex, center_lon: float, center_lat: float) -> tuple:
    """
    Create empty spatial and temporal aggregate datasets for regions without DH areas.
    
    Parameters
    ----------
    snapshots : pd.DatetimeIndex
        Time snapshots for the spatial aggregate
    center_lon : float
        Longitude of region center
    center_lat : float
        Latitude of region center
        
    Returns
    -------
    tuple
        Tuple of (spatial_aggregate, temporal_aggregate) datasets
    """
    spatial_aggregate = xr.Dataset(
        data_vars={
            "total_power": xr.DataArray(
                np.zeros(len(snapshots)),
                dims=["time"],
                coords={"time": snapshots},
            ),
            "average_temperature": xr.DataArray(
                np.zeros(len(snapshots)),
                dims=["time"],
                coords={"time": snapshots},
            ),
        }
    )

    temporal_aggregate = xr.Dataset(
        data_vars={
            "total_energy": xr.DataArray(
                [[0.0]],
                dims=["longitude", "latitude"],
                coords={"longitude": [center_lon], "latitude": [center_lat]},
            ),
            "average_temperature": xr.DataArray(
                [[0.0]],
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
    if region.geometry.is_empty.any():
        # Get the center of the original region (before intersection)
        region_center = original_region.to_crs("EPSG:4326").centroid.iloc[0]
        center_lon = region_center.x
        center_lat = region_center.y
        
        spatial_aggregate, temporal_aggregate = _create_empty_datasets(
            snapshots, center_lon, center_lat
        )
        
        return {
            "spatial aggregate": spatial_aggregate,
            "temporal aggregate": temporal_aggregate,
        }

    # Process region with valid DH area intersection
    # Get bounding box for region
    minx, miny, maxx, maxy = region.total_bounds

    logging.info("Loading river-water data and boxing to region bounds...")

    # Load and preprocess river discharge data
    river_discharge = (
        xr.open_dataset(
            river_discharge_fn,
            chunks={"time": "auto", "lat": "auto", "lon": "auto"},
        )["dis"]
        .rename({"lat": "latitude", "lon": "longitude"})
        .rio.write_crs("EPSG:4326")
        .rio.clip_box(minx, miny, maxx, maxy)
        .rio.reproject("EPSG:3035")
    )

    # Load and preprocess ambient temperature data
    ambient_temperature = (
        xr.open_dataset(
            ambient_temperature_fn,
            chunks={"time": "auto", "lat": "auto", "lon": "auto"},
        )["ta6"]
        .rename({"lat": "latitude", "lon": "longitude"})
        .rio.write_crs("EPSG:4326")
        .rio.clip_box(minx, miny, maxx, maxy)
        .rio.reproject("EPSG:3035")
    )

    # Use projected CRS for spatial calculations
    region = region.to_crs("EPSG:3035")

    logging.info("Approximating river-heat potential...")
    # Instantiate and run the river water heat approximator
    river_water_heat_approximator = RiverWaterHeatApproximator(
        volume_flow=river_discharge,
        ambient_temperature=ambient_temperature,
        region=region,
    )
    spatial_aggregate = river_water_heat_approximator.get_spatial_aggregate().compute()
    temporal_aggregate = (
        river_water_heat_approximator.get_temporal_aggregate()
        .rio.reproject("EPSG:4326")
        .rename({"x": "longitude", "y": "latitude"})
        .compute()
    )

    return {
        "spatial aggregate": spatial_aggregate,
        # temporal aggregate is only used for plotting/analysis
        "temporal aggregate": temporal_aggregate,
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


    # Set up Dask cluster for parallel computation
    cluster = LocalCluster(
        n_workers=int(snakemake.threads),
        threads_per_worker=1,
        memory_limit=f"{snakemake.resources.mem_mb / snakemake.threads}MB",
    )
    client = Client(cluster)

    # Process each region in parallel
    futures = []
    for region_name in regions_onshore.index:
        logging.info(f"Processing region {region_name}")
        region = gpd.GeoSeries(regions_onshore.loc[region_name].copy(deep=True))
        futures.append(
            get_regional_result(
                river_discharge_fn=snakemake.input.hera_river_discharge,
                ambient_temperature_fn=snakemake.input.hera_ambient_temperature,
                region=region,
                dh_areas=dh_areas,
                snapshots=snapshots,
            ),
        )

    # Gather results from all regions
    results = client.gather(futures)

    # Build DataFrame of total power for each region
    power = pd.DataFrame(
        {
            region_name: res["spatial aggregate"]["total_power"].to_pandas()
            for region_name, res in zip(regions_onshore.index, results)
        }
    ).dropna()

    # Align power DataFrame to simulation snapshots
    power = power.reindex(snapshots, method="nearest")
    power.to_csv(snakemake.output.heat_source_power)

    # Concatenate average temperature for all regions
    temperature = xr.concat(
        [res["spatial aggregate"]["average_temperature"] for res in results], dim="name"
    ).assign_coords(name=regions_onshore.index).dropna(dim="time")

    # Align temperature to simulation snapshots
    temperature = temperature.sel(time=snapshots, method="nearest").assign_coords(
        time=snapshots
    )
    temperature.to_netcdf(snakemake.output.heat_source_temperature)

    # Merge the temporal aggregate results for energy and temperature
    xr.concat(
        [res["temporal aggregate"]["total_energy"] for res in results],
        dim=regions_onshore.index,
    ).to_netcdf(snakemake.output.heat_source_energy_temporal_aggregate)

    xr.concat(
        [res["temporal aggregate"]["average_temperature"] for res in results],
        dim=regions_onshore.index,
    ).to_netcdf(snakemake.output.heat_source_temperature_temporal_aggregate)
