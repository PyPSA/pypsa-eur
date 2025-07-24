# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import logging

import geopandas as gpd
import pandas as pd
import xarray as xr
import numpy as np
from _helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
    update_config_from_wildcards,
)
from approximators.sea_water_heat_approximator import SeaWaterHeatApproximator
from dask.distributed import Client, LocalCluster

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
            "average_temperature": xr.DataArray(
                np.zeros(len(snapshots)),
                dims=["time"],
                coords={"time": snapshots},
            ),
        }
    )

    temporal_aggregate = xr.Dataset(
        data_vars={
            "average_temperature": xr.DataArray(
                [[0.0]],
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
    # Get bounds for initial clip_box (more efficient)
    minx, miny, maxx, maxy = region.total_bounds

    water_temperature = (
        xr.open_dataset(
            seawater_temperature_fn,
            chunks={
                "time": "auto",
                "latitude": "auto",
                "longitude": "auto",
                "depth": 1,
            },
            # engine="rasterio",
        )["thetao"]
        .mean(dim="depth")
        .rio.write_crs("EPSG:4326")
        .rio.clip_box(minx, miny, maxx, maxy)
        .rio.reproject("EPSG:3035")
    )

    region = region.to_crs("EPSG:3035")

    seawater_heat_approximator = SeaWaterHeatApproximator(
        water_temperature=water_temperature,
        region=region,
    )

    return {
        "spatial aggregate": seawater_heat_approximator.get_spatial_aggregate().compute(),
        # temporal aggregate is only used for plotting/analysis
        "temporal aggregate": seawater_heat_approximator.get_temporal_aggregate()
        .rio.reproject("EPSG:4326")
        .rename({"x": "longitude", "y": "latitude"})
        .compute(),
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

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    snapshots: pd.DatetimeIndex = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    regions_onshore = gpd.read_file(snakemake.input["regions_onshore"])
    regions_onshore.set_index("name", inplace=True)
    regions_onshore.set_crs("EPSG:4326", inplace=True)

    dh_areas = gpd.read_file(snakemake.input["dh_areas"]).to_crs("EPSG:3035")
    dh_areas["geometry"] = dh_areas.geometry.buffer(snakemake.params.dh_area_buffer)
    dh_areas = dh_areas.to_crs("EPSG:4326")

    cluster = LocalCluster(
        n_workers=int(snakemake.threads),
        threads_per_worker=1,
        memory_limit=f"{snakemake.resources.mem_mb / snakemake.threads}MB",
    )
    client = Client(cluster)

    futures = []
    for region_name in regions_onshore.index:
        logging.info(f"Processing region {region_name}")
        region = gpd.GeoSeries(regions_onshore.loc[region_name].copy(deep=True))
        futures.append(
            get_regional_result(
                seawater_temperature_fn=snakemake.input["seawater_temperature"],
                region=region,
                dh_areas=dh_areas,
                snapshots=snapshots,
            )
        )

    results = client.gather(futures)

    temperature = xr.concat(
        [res["spatial aggregate"]["average_temperature"] for res in results], dim="name"
    ).assign_coords(name=regions_onshore.index).dropna(dim="time")

    temperature = temperature.sel(time=snapshots, method="nearest").assign_coords(
        time=snapshots
    )
    temperature.to_netcdf(snakemake.output.heat_source_temperature)

    # Merge the temporal aggregate results
    xr.concat(
        [res["temporal aggregate"]["average_temperature"] for res in results],
        dim=regions_onshore.index,
    ).to_netcdf(snakemake.output.heat_source_temperature_temporal_aggregate)
