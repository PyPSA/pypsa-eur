# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import logging

import geopandas as gpd
import pandas as pd
import xarray as xr
from _helpers import (
    configure_logging,
    get_snapshots,
    set_scenario_config,
    update_config_from_wildcards,
)
from dask.distributed import Client, LocalCluster

logger = logging.getLogger(__name__)


from approximators.river_water_heat_approximator import RiverWaterHeatApproximator


def get_regional_result(
    river_discharge_fn: xr.DataArray,
    ambient_temperature_fn: xr.DataArray,
    region: gpd.GeoSeries,
    dh_areas: gpd.GeoDataFrame,
) -> dict:
    # Clip the region to the district heating areas
    region.geometry = gpd.overlay(
        region.to_frame(),
        dh_areas,
        how="intersection",
    ).union_all()
    # Get bounds for initial clip_box
    minx, miny, maxx, maxy = region.total_bounds

    logging.info("Loading river-water data and boxing to region bounds...")

    river_discharge = (
        xr.open_dataset(
            river_discharge_fn,
            chunks={"time": "auto", "lat": "auto", "lon": "auto"},
            # engine="rasterio",
        )["dis"]
        .rename({"lat": "latitude", "lon": "longitude"})
        .rio.write_crs("EPSG:4326")
        .rio.clip_box(minx, miny, maxx, maxy)
        .rio.reproject("EPSG:3035")
    )

    ambient_temperature = (
        xr.open_dataset(
            ambient_temperature_fn,
            chunks={"time": "auto", "lat": "auto", "lon": "auto"},
            # engine="rasterio",
        )["ta6"]
        .rename({"lat": "latitude", "lon": "longitude"})
        .rio.write_crs("EPSG:4326")
        .rio.clip_box(minx, miny, maxx, maxy)
        .rio.reproject("EPSG:3035")
    )

    # Use EPSG:3035 for all calculations to use buffers/distances etc. in meters
    region = region.to_crs("EPSG:3035")

    logging.info("Approximating river-heat potential...")
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

    configure_logging(snakemake)
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    snapshots: pd.DatetimeIndex = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )

    regions_onshore = gpd.read_file(snakemake.input["regions_onshore"])
    regions_onshore.set_index("name", inplace=True)
    regions_onshore = regions_onshore.to_crs("EPSG:4326")

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
            # results.append(
            get_regional_result(
                river_discharge_fn=snakemake.input.hera_river_discharge,
                ambient_temperature_fn=snakemake.input.hera_ambient_temperature,
                region=region,
                dh_areas=dh_areas,
            ),
        )

    results = client.gather(futures)

    power = pd.DataFrame(
        {
            region_name: res["spatial aggregate"]["total_power"].to_pandas()
            for res in results
            for region_name, res in zip(regions_onshore.index, results)
        }
    )

    power = power.reindex(snapshots, method="nearest")
    power.to_csv(snakemake.output.heat_source_power)

    temperature = xr.concat(
        [res["spatial aggregate"]["average_temperature"] for res in results], dim="name"
    ).assign_coords(name=regions_onshore.index)

    temperature = temperature.sel(time=snapshots, method="nearest").assign_coords(
        time=snapshots
    )
    temperature.to_netcdf(snakemake.output.heat_source_temperature)

    # Merge the temporal aggregate results
    xr.concat(
        [res["temporal aggregate"]["total_energy"] for res in results],
        dim=regions_onshore.index,
    ).to_netcdf(snakemake.output.heat_source_energy_temporal_aggregate)

    xr.concat(
        [res["temporal aggregate"]["average_temperature"] for res in results],
        dim=regions_onshore.index,
    ).to_netcdf(snakemake.output.heat_source_temperature_temporal_aggregate)
