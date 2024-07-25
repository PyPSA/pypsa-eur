# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import numpy as np
import xarray as xr
from _helpers import set_scenario_config
from CentralHeatingCopApproximator import CentralHeatingCopApproximator
from DecentralHeatingCopApproximator import DecentralHeatingCopApproximator

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
        source_inlet_temperature_celsius = xr.open_dataarray(
            snakemake.input[f"temp_{source_type}_total"]
        )

        # Approximate COP for decentral (individual) heating
        cop_individual_heating = DecentralHeatingCopApproximator(
            forward_temperature_celsius=snakemake.params.heat_pump_sink_T_decentral_heating,
            source_inlet_temperature_celsius=source_inlet_temperature_celsius,
            source_type=source_type,
        ).approximate_cop()
        cop_individual_heating.to_netcdf(
            snakemake.output[f"cop_{source_type}_decentral_heating"]
        )

        # Approximate COP for central (district) heating
        cop_central_heating = CentralHeatingCopApproximator(
            forward_temperature_celsius=snakemake.params.forward_temperature_central_heating,
            return_temperature_celsius=snakemake.params.return_temperature_central_heating,
            source_inlet_temperature_celsius=source_inlet_temperature_celsius,
            source_outlet_temperature_celsius=source_inlet_temperature_celsius
            - snakemake.params.heat_source_cooling_central_heating,
        ).approximate_cop()
        cop_central_heating.to_netcdf(
            snakemake.output[f"cop_{source_type}_central_heating"]
        )
