# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Approximate heat pump coefficient-of-performance (COP) profiles for different
heat sources and systems.

For central heating, this is based on Jensen et al. (2018) (c.f. `CentralHeatingCopApproximator <CentralHeatingCopApproximator.py>`_) and for decentral heating, the approximation is based on Staffell et al. (2012) (c.f. `DecentralHeatingCopApproximator <DecentralHeatingCopApproximator.py>`_).

Relevant Settings
-----------------

.. code:: yaml
    sector:
        heat_pump_sink_T_decentral_heating:
        district_heating:
            forward_temperature:
            return_temperature:
            heat_source_cooling:
            heat_pump_cop_approximation:
                refrigerant:
                heat_exchanger_pinch_point_temperature_difference
                isentropic_compressor_efficiency:
                heat_loss:
            heat_pump_sources:
                urban central:
                urban decentral:
                rural:
    snapshots:

Inputs
------
- `resources/<run_name>/regions_onshore.geojson`: Onshore regions
- `resources/<run_name>/temp_soil_total`: Ground temperature
- `resources/<run_name>/temp_air_total`: Air temperature

Outputs
-------
- `resources/<run_name>/cop_profiles.nc`: Heat pump coefficient-of-performance (COP) profiles
"""

import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import set_scenario_config
from CentralHeatingCopApproximator import CentralHeatingCopApproximator
from DecentralHeatingCopApproximator import DecentralHeatingCopApproximator

from scripts.definitions.heat_system_type import HeatSystemType


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


def get_cop(
    heat_system_type: str,
    heat_source: str,
    source_inlet_temperature_celsius: xr.DataArray,
    forward_temperature_by_node_and_time: xr.DataArray = None,
    return_temperature_by_node_and_time: xr.DataArray = None,
) -> xr.DataArray:
    """
    Calculate the coefficient of performance (COP) for a heating system.

    Parameters
    ----------
    heat_system_type : str
        The type of heating system.
    heat_source : str
        The heat source used in the heating system.
    source_inlet_temperature_celsius : xr.DataArray
        The inlet temperature of the heat source in Celsius.

    Returns
    -------
    xr.DataArray
        The calculated coefficient of performance (COP) for the heating system.
    """
    if HeatSystemType(heat_system_type).is_central:
        return CentralHeatingCopApproximator(
            forward_temperature_celsius=forward_temperature_by_node_and_time,
            return_temperature_celsius=return_temperature_by_node_and_time,
            source_inlet_temperature_celsius=source_inlet_temperature_celsius,
            source_outlet_temperature_celsius=source_inlet_temperature_celsius
            - snakemake.params.heat_source_cooling_central_heating,
        ).approximate_cop()

    else:
        return DecentralHeatingCopApproximator(
            forward_temperature_celsius=snakemake.params.heat_pump_sink_T_decentral_heating,
            source_inlet_temperature_celsius=source_inlet_temperature_celsius,
            source_type=heat_source,
        ).approximate_cop()


def get_country_from_node_name(node_name: str) -> str:
    return node_name[:2]


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
    forward_temperature_central_heating_by_node_and_time: xr.DataArray = (
        map_temperature_dict_to_onshore_regions(
            supply_temperature_by_country=snakemake.params.forward_temperature_central_heating,
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
    cop_all_system_types = []
    for heat_system_type, heat_sources in snakemake.params.heat_pump_sources.items():
        cop_this_system_type = []
        for heat_source in heat_sources:
            source_inlet_temperature_celsius = xr.open_dataarray(
                snakemake.input[f"temp_{heat_source.replace('ground', 'soil')}_total"]
            )
            cop_da = get_cop(
                heat_system_type=heat_system_type,
                heat_source=heat_source,
                source_inlet_temperature_celsius=source_inlet_temperature_celsius,
                forward_temperature_by_node_and_time=forward_temperature_central_heating_by_node_and_time,
                return_temperature_by_node_and_time=return_temperature_central_heating_by_node_and_time,
            )
            cop_this_system_type.append(cop_da)
        cop_all_system_types.append(
            xr.concat(
                cop_this_system_type, dim=pd.Index(heat_sources, name="heat_source")
            )
        )

    cop_dataarray = xr.concat(
        cop_all_system_types,
        dim=pd.Index(snakemake.params.heat_pump_sources.keys(), name="heat_system"),
    )

    cop_dataarray.to_netcdf(snakemake.output.cop_profiles)
