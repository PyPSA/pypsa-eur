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
            clusters=48,
        )

    set_scenario_config(snakemake)

    # map forward and return temperatures specified on country-level to onshore regions
    regions_onshore = gpd.read_file(snakemake.input.regions_onshore)["name"]
    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots)
    central_heating_forward_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    central_heating_return_temperature: xr.DataArray = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
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
                forward_temperature_by_node_and_time=central_heating_forward_temperature,
                return_temperature_by_node_and_time=central_heating_return_temperature,
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
