# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Generate PTES temperature profiles.

This script calculates forward and return temperature adjustments and models
the base temperature profile for Pit Thermal Energy Storage (PTES).

Relevant Settings
-----------------
.. code:: yaml
    storage:
        PTES:
            max_temperature: 90 # Stimmt nicht
            min_temperature: 58 # Stimmt nicht, die min_temperatur wird durch die return_temp festgelegt

Inputs
------
- `resources/<run_name>/forward_temperature.nc`: Forward temperature profile
- `resources/<run_name>/return_temperature.nc`: Return temperature profile

Outputs
-------
- `resources/<run_name>/simplified_temperature_model.nc`:
"""

# beachte pinch point ist 5K in CentralHeatingCopApproximator

import numpy as np
import pandas as pd
import xarray as xr
from build_TES_temperature_profiles import PTESTemperatureApproximator

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from _helpers import set_scenario_config

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tes_storage_temperature_profile",
            clusters=6,
            planning_horizons="2030",
        )

    set_scenario_config(snakemake)

    # Load inputs
    forward_temp = xr.open_dataarray(snakemake.input.central_heating_forward_temperature_profiles)
    return_temp = xr.open_dataarray(snakemake.input.central_heating_return_temperature_profiles)

    # Load the snapshots
    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots)

    # Load configuration
    max_PTES_temperature = snakemake.params.max_PTES_temperature

    # Initialize the PTES Temperature Approximator
    ptes_approximator = PTESTemperatureApproximator(
        forward_temperature_celsius=forward_temp,
        return_temperature_celsius=return_temp,
        max_PTES_temperature=max_PTES_temperature,
        snapshots=snapshots,
    )

    # Calculate simplified temperature model
    ltes_top_layer_temperature_profiles = ptes_approximator.simplified_top_layer_temperature_model()
    ltes_bottom_layer_temperature_profiles = ptes_approximator.ptes_bottom_temperature()

    # Save outputs
    ltes_top_layer_temperature_profiles.to_netcdf(
        snakemake.output.ltes_top_layer_temperature_profiles
    )
    ltes_bottom_layer_temperature_profiles.to_netcdf(
        snakemake.output.ltes_bottom_layer_temperature_profiles
    )
