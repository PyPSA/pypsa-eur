import geopandas as gpd
import pandas as pd
from dask.distributed import Client, LocalCluster
import shapely
import xarray as xr

import logging
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
    get_snapshots,
)

logger = logging.getLogger(__name__)


from seawater_thermal_approximator import SeawaterThermalApproximator


def get_regional_result(
    seawater_data_fn: str, geometry: shapely.geometry.polygon.Polygon
) -> dict:

    ds = xr.open_dataset(
        seawater_data_fn,
        chunks={"time": 8760, "latitude": 50, "longitude": 50},
        decode_coords=["time", "latitude", "longitude"],
        mode="r",
    ).sel(
        longitude=slice(geometry.bounds[0], geometry.bounds[2]),
        latitude=slice(
            geometry.bounds[1],
            geometry.bounds[3],
        ),
    )

    return SeawaterThermalApproximator(data=ds, geometry=geometry).results.compute()


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

    noprogress = snakemake.config["run"].get("disable_progressbar", True)
    noprogress = noprogress or not snakemake.config["atlite"]["show_progress"]

    regions_onshore = gpd.read_file(snakemake.input["regions_onshore"])
    regions_onshore.set_index("name", inplace=True)
    regions_onshore["name"] = regions_onshore.index

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
                seawater_data_fn=snakemake.input["seawater_data"], geometry=geometry
            )
        )

    results = client.gather(futures)

    power = pd.DataFrame(
        {
            region_name: res["total_power"].to_pandas()
            for res in results
            for region_name, res in zip(regions_onshore.index, results)
        }
    )

    power = power.reindex(snapshots, method="nearest")
    power.to_csv(snakemake.output.heat_source_power)

    temperature = xr.concat(
        [res["average_temperature"] for res in results], dim="name"
    ).assign_coords(name=regions_onshore.index)

    temperature = temperature.sel(time=snapshots, method="nearest").assign_coords(time=snapshots)
    temperature.to_netcdf(snakemake.output.heat_source_temperature)
