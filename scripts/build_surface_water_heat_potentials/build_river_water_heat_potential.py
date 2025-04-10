# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
import logging

import geopandas as gpd
import pandas as pd
import shapely
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
    river_discharge: xr.DataArray,
    ambient_temperature: xr.DataArray,
    geometry: shapely.geometry.polygon.Polygon,
) -> dict:
    river_water_heat_approximator = RiverWaterHeatApproximator(
        volume_flow=river_discharge,
        ambient_temperature=ambient_temperature,
        region_geometry=geometry,
    )
    spatial_aggregate = river_water_heat_approximator.get_spatial_aggregate().compute()
    temporal_aggregate = (
        river_water_heat_approximator.get_temporal_aggregate().compute()
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
    regions_onshore["name"] = regions_onshore.index

    cluster = LocalCluster(
        n_workers=int(snakemake.threads),
        threads_per_worker=1,
        memory_limit=f"{snakemake.resources.mem_mb / snakemake.threads}MB",
    )
    client = Client(cluster)

    river_discharge = (
        xr.open_dataset(
            snakemake.input.hera_river_discharge,
            chunks={"time": -1, "lat": 100, "lon": 100},
            decode_coords=["time", "lat", "lon"],
            mode="r",
        )["dis"]
        .sortby(["time", "lat", "lon"])
        .rename({"lon": "longitude", "lat": "latitude"})
    )

    ambient_temperature = (
        xr.open_dataset(
            snakemake.input.hera_ambient_temperature,
            chunks={"time": -1, "lat": 100, "lon": 100},
            decode_coords=["time", "lat", "lon"],
            mode="r",
        )["ta6"]
        .sortby(["time", "lat", "lon"])
        .rename({"lon": "longitude", "lat": "latitude"})
    )

    futures = []
    for region_name in regions_onshore.index:
        geometry = regions_onshore.loc[region_name].geometry
        futures.append(
            get_regional_result(
                river_discharge=river_discharge.sel(
                    longitude=slice(geometry.bounds[0], geometry.bounds[2]),
                    latitude=slice(
                        geometry.bounds[1],
                        geometry.bounds[3],
                    ),
                ),
                ambient_temperature=ambient_temperature.sel(
                    longitude=slice(geometry.bounds[0], geometry.bounds[2]),
                    latitude=slice(
                        geometry.bounds[1],
                        geometry.bounds[3],
                    ),
                ),
                geometry=geometry,
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
