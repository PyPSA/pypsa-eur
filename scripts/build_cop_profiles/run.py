# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import numpy as np
import xarray as xr
from _helpers import set_scenario_config, get_country_from_node_name
from CentralHeatingCopApproximator import CentralHeatingCopApproximator
from DecentralHeatingCopApproximator import DecentralHeatingCopApproximator


def map_temperature_dict_to_onshore_regions(
    temperature_dict: dict, onshore_regions: xr.DataArray, snapshots: xr.DataArray
) -> xr.DataArray:
    """
    Map dictionary of temperatures to onshore regions.

    Parameters:
    ----------
    temperature_dict : dictionary
        Dictionary with temperatures as values and country keys as keys. One key must be named "default"
    onshore_regions : xr.DataArray
        Names of onshore regions

    Returns:
    -------
    xr.DataArray
        The dictionary values mapped to onshore regions with onshore regions as coordinates.
    """
    return xr.DataArray(
        [
            [
                (
                    temperature_dict[get_country_from_node_name(node_name)]
                    if get_country_from_node_name(node_name) in temperature_dict.keys()
                    else temperature_dict["default"]
                )
                for node_name in onshore_regions["name"].values
            ]
            # pass both nodes and snapshots as dimensions to preserve correct data structure
            for _ in snapshots["time"].values
        ],
        dims=["time", "name"],
        coords={"time": snapshots["time"], "name": onshore_regions["name"]},
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

    for source_type in ["air", "soil"]:
        # source inlet temperature (air/soil) is based on weather data
        source_inlet_temperature_celsius: xr.DataArray = xr.open_dataarray(
            snakemake.input[f"temp_{source_type}_total"]
        )

        # Approximate COP for decentral (individual) heating
        cop_individual_heating: xr.DataArray = DecentralHeatingCopApproximator(
            forward_temperature_celsius=snakemake.params.heat_pump_sink_T_decentral_heating,
            source_inlet_temperature_celsius=source_inlet_temperature_celsius,
            source_type=source_type,
        ).approximate_cop()
        cop_individual_heating.to_netcdf(
            snakemake.output[f"cop_{source_type}_decentral_heating"]
        )

        # map forward and return temperatures specified on country-level to onshore regions
        onshore_regions: xr.DataArray = source_inlet_temperature_celsius["name"]
        forward_temperature_central_heating: xr.DataArray = (
            map_temperature_dict_to_onshore_regions(
                temperature_dict=snakemake.params.forward_temperature_central_heating,
                onshore_regions=onshore_regions,
                snapshots=source_inlet_temperature_celsius["time"]
            )
        )
        return_temperature_central_heating: xr.DataArray = (
            map_temperature_dict_to_onshore_regions(
                temperature_dict=snakemake.params.return_temperature_central_heating,
                onshore_regions=onshore_regions,
                snapshots=source_inlet_temperature_celsius["time"]
            )
        )

        # Approximate COP for central (district) heating
        cop_central_heating: xr.DataArray = CentralHeatingCopApproximator(
            forward_temperature_celsius=forward_temperature_central_heating,
            return_temperature_celsius=return_temperature_central_heating,
            source_inlet_temperature_celsius=source_inlet_temperature_celsius,
            source_outlet_temperature_celsius=source_inlet_temperature_celsius
            - snakemake.params.heat_source_cooling_central_heating,
        ).approximate_cop()
        cop_central_heating.to_netcdf(
            snakemake.output[f"cop_{source_type}_central_heating"]
        )
