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


from approximators.sea_water_heat_approximator import SeaWaterHeatApproximator


def get_regional_result(
    seawater_temperature_fn: str, geometry: shapely.geometry.polygon.Polygon
) -> dict:
    ds = xr.open_dataset(
        seawater_temperature_fn,
        chunks={"time": 8760, "latitude": 50, "longitude": 50},
        decode_coords=["time", "latitude", "longitude"],
        mode="r",
    ).sortby(["time", "latitude", "longitude"]).sel(
        longitude=slice(geometry.bounds[0], geometry.bounds[2]),
        latitude=slice(
            geometry.bounds[1],
            geometry.bounds[3],
        ),
    )

    seawater_heat_approximator = SeaWaterHeatApproximator(
        water_temperature=ds["thetao"].mean(dim="depth"),
        region_geometry=geometry,
    )

    return {
        "spatial aggregate": seawater_heat_approximator.get_spatial_aggregate().compute(),
        # temporal aggregate is only used for plotting/analysis
        "temporal aggregate": seawater_heat_approximator.get_temporal_aggregate().compute(),
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

    cluster = LocalCluster(
        n_workers=int(snakemake.threads),
        threads_per_worker=1,
        memory_limit=f"{snakemake.resources.mem_mb / snakemake.threads}MB",
    )
    client = Client(cluster)

    futures = []
    for region_name in regions_onshore.index:
        geometry = regions_onshore.loc[region_name].geometry
        futures.append(
            get_regional_result(
                seawater_temperature_fn=snakemake.input["seawater_temperature"],
                geometry=geometry,
            )
        )

    results = client.gather(futures)

    temperature = xr.concat(
        [res["spatial aggregate"]["average_temperature"] for res in results], dim="name"
    ).assign_coords(name=regions_onshore.index)

    temperature = temperature.sel(time=snapshots, method="nearest").assign_coords(
        time=snapshots
    )
    temperature.to_netcdf(snakemake.output.heat_source_temperature)

    # Merge the temporal aggregate results
    xr.concat(
        [res["temporal aggregate"]["average_temperature"] for res in results],
        dim=regions_onshore.index,
    ).to_netcdf(snakemake.output.heat_source_temperature_temporal_aggregate)
