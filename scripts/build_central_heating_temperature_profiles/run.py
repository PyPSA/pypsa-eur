# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Approximate district heating forward and return temperature profiles based on
ambient temperature. The method is based on a reference curve from Pieper et
al. 2019, where for ambient temperatures below 0C, the highest possible forward
temperature is assumed and vice versa for temperatures above 10C. Between these
threshold levels, forward temperatures are linearly interpolated.

By default, temperature levels are increased for non-Scandinavian countries.
The default ratios between min. and max. forward temperatures is based on AGFW-Hauptbericht 2022.

Relevant Settings
-----------------

.. code:: yaml
    sector:
        district_heating:
            max_forward_temperature:
            min_forward_temperature:
            return_temperature:
Inputs
------
- `resources/<run_name>/temp_air_total`: Air temperature

Outputs
-------
- `resources/<run_name>/central_heating_temperature_profiles.nc`:

References
----------
- Pieper, et al. (2019): "Assessment of a combination of three heat sources for heat pumps to supply district heating" (https://doi.org/10.1016/j.energy.2019.03.165).
- AGFW (2022): "Hauptbericht 2022" (https://www.agfw.de/zahlen-und-statistiken/agfw-hauptbericht)
"""

import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import set_scenario_config
from central_heating_temperature_approximator import (
    CentralHeatingTemperatureApproximator,
)


def get_country_from_node_name(node_name: str) -> str:
    return node_name[:2]


def map_temperature_dict_to_onshore_regions(
    supply_temperature_by_country: dict,
    regions_onshore: pd.Index,
    snapshots: pd.DatetimeIndex,
) -> xr.DataArray:
    """
    Map dictionary of temperatures to onshore regions.

    Parameters:
    ----------
    supply_temperature_by_country : dictionary
        Dictionary with temperatures as values and country keys as keys. One key must be named "default"
    regions_onshore : pd.Index
        Names of onshore regions
    snapshots : pd.DatetimeIndex
        Time stamps

    Returns:
    -------
    xr.DataArray
        The dictionary values mapped to onshore regions with onshore regions as coordinates.
    """
    return xr.DataArray(
        [
            [
                (
                    supply_temperature_by_country[get_country_from_node_name(node_name)]
                    if get_country_from_node_name(node_name)
                    in supply_temperature_by_country.keys()
                    else supply_temperature_by_country["default"]
                )
                for node_name in regions_onshore.values
            ]
            # pass both nodes and snapshots as dimensions to preserve correct data structure
            for _ in snapshots
        ],
        dims=["time", "name"],
        coords={"time": snapshots, "name": regions_onshore},
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            simpl="",
            clusters=48,
        )

    set_scenario_config(snakemake)

    # map forward and return temperatures specified on country-level to onshore regions
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore)["name"]
    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots)
    max_forward_temperature_central_heating_by_node_and_time: xr.DataArray = (
        map_temperature_dict_to_onshore_regions(
            supply_temperature_by_country=snakemake.params.max_forward_temperature_central_heating,
            regions_onshore=regions_onshore,
            snapshots=snapshots,
        )
    )
    min_forward_temperature_central_heating_by_node_and_time: xr.DataArray = (
        map_temperature_dict_to_onshore_regions(
            supply_temperature_by_country=snakemake.params.min_forward_temperature_central_heating,
            regions_onshore=regions_onshore,
            snapshots=snapshots,
        )
    )
    return_temperature_central_heating_by_node_and_time: xr.DataArray = (
        map_temperature_dict_to_onshore_regions(
            supply_temperature_by_country=snakemake.params.return_temperature_central_heating,
            regions_onshore=regions_onshore,
            snapshots=snapshots,
        )
    )

    central_heating_temperature_approximator = CentralHeatingTemperatureApproximator(
        ambient_temperature=xr.open_dataarray(snakemake.input.temp_air_total),
        max_forward_temperature=max_forward_temperature_central_heating_by_node_and_time,
        min_forward_temperature=min_forward_temperature_central_heating_by_node_and_time,
        fixed_return_temperature=return_temperature_central_heating_by_node_and_time,
        lower_threshold_ambient_temperature=snakemake.params.lower_threshold_ambient_temperature,
        upper_threshold_ambient_temperature=snakemake.params.upper_threshold_ambient_temperature,
        rolling_window_ambient_temperature=snakemake.params.rolling_window_ambient_temperature,
    )

    central_heating_temperature_approximator.forward_temperature.to_netcdf(
        snakemake.output.central_heating_forward_temperature_profiles
    )
    central_heating_temperature_approximator.return_temperature.to_netcdf(
        snakemake.output.central_heating_return_temperature_profiles
    )
