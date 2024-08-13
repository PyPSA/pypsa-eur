# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import sys

import numpy as np
import pandas as pd
import xarray as xr
from _helpers import set_scenario_config
from CentralHeatingCopApproximator import CentralHeatingCopApproximator
from DecentralHeatingCopApproximator import DecentralHeatingCopApproximator

from scripts.definitions.heat_system_type import HeatSystemType

sys.path.append("..")


def get_cop(
    heat_system_type: str,
    heat_source: str,
    source_inlet_temperature_celsius: xr.DataArray,
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
            forward_temperature_celsius=snakemake.params.forward_temperature_central_heating,
            return_temperature_celsius=snakemake.params.return_temperature_central_heating,
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


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            simpl="",
            clusters=48,
        )

    set_scenario_config(snakemake)

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
