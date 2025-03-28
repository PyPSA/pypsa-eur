# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build heat source power and temperature time series from river water. The outputs are indexed by `snapshots` and `name` (the name of the onshore region).
The heat source power is in MW and the heat source temperature is in Kelvin.

River power and temperature are approximated using the implementation in `RiverThermalApproximator <RiverThermalApproximator.py>`_.

There, river water temperature is approximated from the smoothed ambient temperature, based on Triebs & Trsatsaronis 2022: Estimating the local renewable potentials for the transformation of district heating systems, Proceedings of ECOS 2022, pp. 479–490: https://orbit.dtu.dk/en/publications/proceedings-of-ecos-2022-the-35th-international-conference-on-eff.

River power is the product of the river water volume flow extracted, the temperature difference between the river and the minimum inlet temperature, the heat capacity and the density of water.

Relevant Settings
-----------------
.. code:: yaml
    snapshots:
    drop_leap_day:

Inputs
------
- `resources/<run_name>/regions_onshore.geojson`
- `resources/<run_name>/hera_river_discharge.nc`
- `resources/<run_name>/hera_ambient_temperature.nc`

Outputs
-------
- `resources/<run_name>/heat_source_power_river_water_base_s_{clusters}.nc`
- `resources/<run_name>/heat_source_temperature_river_water_base_s_{clusters}.nc`
"""

import pandas as pd
import geopandas as gpd
import xarray as xr
from typing import Tuple
import os
from dask.distributed import Client

from _helpers import set_scenario_config, get_snapshots
from river_thermal_approximator import RiverThermalApproximator

# DEBUG_MODE = True


def chunk_hera_data(
    hera_data: xr.Dataset, chunks: dict = {"time": -1, "lat": "auto", "lon": "auto"}
) -> xr.Dataset:
    """
    Chunk the hera data.

    Parameters
    ----------
    hera_data : xr.DataSet
        The data from hera.
    chunks : dict
        The chunks to use. Default is {"time": -1, "lat": "auto", "lon": "auto"}.

    Returns
    -------
    xr.Dataset
        The chunked data.
    """
    return hera_data.chunk(chunks)


def box_hera_data_to_onshore_regions(
    hera_data: xr.Dataset, regions_onshore: gpd.GeoDataFrame
) -> xr.Dataset:
    """
    Get the data from hera in the bounding box of the onshore region.

    Parameters
    ----------
    hera_data : xr.DataSet
        The data from hera.
    regions_onshore : gpd.GeoDataFrame
        The onshore regions.

    Returns
    -------
    xr.Dataset
        The data from hera in the bounding box of the onshore region.
    """

    hera_data = hera_data.sortby(["time", "lat", "lon"])

    return hera_data.sel(
        lat=slice(regions_onshore.bounds.miny.min(), regions_onshore.bounds.maxy.max()),
        lon=slice(regions_onshore.bounds.minx.min(), regions_onshore.bounds.maxx.max()),
    )


def run_approximator(
        geometry: gpd.GeoDataFrame,
        hera_river_discharge: xr.Dataset,
        hera_ambient_temperature: xr.Dataset,
        snapshots: pd.DatetimeIndex,
):

    # build approximator
    approximator = RiverThermalApproximator(
            geometry=geometry,
            hera_river_discharge=hera_river_discharge,
            hera_ambient_temperature=hera_ambient_temperature,
            target_snapshots=snapshots,
        )
    
    # run approximator (caches total_river_power)
    _ = approximator.total_river_power

    return approximator

def get_heat_source_power_and_temperature(
    snapshots: pd.DatetimeIndex,
    regions_onshore: gpd.GeoDataFrame,
    hera_river_discharge: xr.Dataset,
    hera_ambient_temperature: xr.Dataset,
    debug_mode: bool = False,
    client: Client = None,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Get the heat source power and temperature per clustered region.

    Parameters
    ----------
    snapshots : pd.DatetimeIndex
        The snapshots.
    regions_onshore : gpd.GeoDataFrame
        The onshore regions.
    hera_river_discharge : xr.Dataset
        The river discharge data. Index by "lon", "lat", "time". Contains variable "dis".
    hera_ambient_temperature : xr.Dataset
        The ambient temperature data. Index by "lon", "lat", "time". Contains variable "ta6".

    Returns
    -------
    Tuple[xr.DataArray, xr.DataArray]
        The heat source power and temperature per clustered region
    """

    kwargs = {
        region_name: {
            "geometry": gpd.GeoDataFrame(region, geometry=[region.geometry]),
            "hera_river_discharge": hera_river_discharge,
            "hera_ambient_temperature": hera_ambient_temperature,
            "snapshots": snapshots}
        for region_name, region in regions_onshore.iterrows()
    }

    if client is not None:
        
        approximator_futures = {
            region_name: client.submit(run_approximator, **kwargs[region_name]
            )
            for region_name, _ in regions_onshore.iterrows()
        }
        # Get results
        approximators = {x: future.result() for x, future in approximator_futures.items()}
    else:
        approximators = {
            region_name: run_approximator(**kwargs[region_name])
            for region_name, _ in regions_onshore.iterrows()
        }
    
    heat_source_power = xr.concat(
        [
            appr.total_river_power.expand_dims(name=[region])
            for region, appr in approximators.items()
        ],
        dim="name",
    )

    heat_source_temperature = xr.concat(
        [
            appr.mean_river_temperature.expand_dims(name=[region])
            for region, appr in approximators.items()
        ],
        dim="name",
    )

    # if debug_mode:
    #     river_power_in_best_locations = xr.concat(
    #         [
    #             appr._river_power_in_best_locations.expand_dims(name=[region])
    #             for region, appr in approximators.items()
    #         ],
    #         dim="name",
    #     )

    #     river_power_in_geometry = xr.concat(
    #         [
    #             appr._river_power_in_geometry.expand_dims(name=[region])
    #             for region, appr in approximators.items()
    #         ],
    #         dim="name",
    #     )
    #     river_temperature_in_geometry = xr.concat(
    #         [
    #             appr._river_power_in_geometry.expand_dims(name=[region])
    #             for region, appr in approximators.items()
    #         ],
    #         dim="name",
    #     )
    #     return (
    #         heat_source_power,
    #         heat_source_temperature,
    #         {
    #             "river_power_in_best_locations": river_power_in_best_locations,
    #             "river_power_in_geometry": river_power_in_geometry,
    #             "river_temperature_in_geometry": river_temperature_in_geometry,
    #         },
    #     )

    # else:
    if True:
        return heat_source_power, heat_source_temperature

def debug_file_path():
    """ get everything in the string up until "resources/"""
    return snakemake.output[0].split("resources/")[0] + "debug/"

if __name__ == "__main__":

    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_heat_source_potentials",
            clusters=48,
        )

    set_scenario_config(snakemake)

    nprocesses = int(snakemake.threads)
    noprogress = snakemake.config["run"].get("disable_progressbar", True)
    noprogress = noprogress or not snakemake.config["atlite"]["show_progress"]



    snapshots: pd.DatetimeIndex = get_snapshots(
        snakemake.params.snapshots, snakemake.params.drop_leap_day
    )
    regions_onshore: gpd.GeoDataFrame = gpd.read_file(snakemake.input.regions_onshore)
    regions_onshore.set_index("name", inplace=True)

    hera_river_discharge = box_hera_data_to_onshore_regions(
        hera_data=chunk_hera_data(
            xr.open_dataset(snakemake.input.hera_river_discharge)
        ),
        regions_onshore=regions_onshore,
    )
    hera_ambient_temperature = box_hera_data_to_onshore_regions(
        hera_data=chunk_hera_data(
            xr.open_dataset(snakemake.input.hera_ambient_temperature)
        ),
        regions_onshore=regions_onshore,
    )

    if nprocesses > 1:
        client = Client(n_workers=nprocesses, threads_per_worker=1)
    else:
        client = None

    # if DEBUG_MODE:
        # heat_source_power, heat_source_temperature, debug_data = get_heat_source_power_and_temperature(
        #     snapshots=snapshots,
        #     regions_onshore=regions_onshore,
        #     hera_river_discharge=hera_river_discharge,
        #     hera_ambient_temperature=hera_ambient_temperature,
        #     debug_mode=True,
        # )

        # if not os.path.exists(debug_file_path()):
        #     os.makedirs(debug_file_path())
        # for key, val in debug_data.items():
        #     val.to_netcdf(debug_file_path() + key + ".nc")


    # else:
    if True:
        import cProfile
        #cProfile.run("""
        heat_source_power, heat_source_temperature = get_heat_source_power_and_temperature(
            snapshots=snapshots,
            regions_onshore=regions_onshore,
            hera_river_discharge=hera_river_discharge,
            hera_ambient_temperature=hera_ambient_temperature,
            client=client,
        )
        #""")

    # convert kW to MW
    (heat_source_power / 1000).to_pandas().T.to_csv(snakemake.output.heat_source_power)
    heat_source_temperature.to_netcdf(snakemake.output.heat_source_temperature)
