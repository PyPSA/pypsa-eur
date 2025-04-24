# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Determine PTES top temperature and additional heating requirements for integration.

This script evaluates the district heating forward temperature profile to:
  1. Approximate the top temperature of the Pit Thermal Energy Storage (PTES) system,
     ensuring that the temperature does not exceed the operational limit.
  2. Determine whether additional (after) heating is needed. A binary indicator is generated:
       - 1: The forward temperature is less than or equal to the PTES maximum; direct usage is possible.
       - 0: The forward temperature exceeds the PTES maximum; additional heating (e.g., via a heat pump) is required.

Relevant Settings
-----------------
.. code:: yaml
    storage:
        PTES:
            max_temperature: 90  # Maximum temperature (°C) the PTES can deliver directly

Inputs
------
- `resources/<run_name>/forward_temperature.nc`
    Forward temperature profile (in Celsius) for the district heating network.

Outputs
-------
- `resources/<run_name>/ptes_top_temperature.nc`
    NetCDF file containing the PTES top temperature profile (clipped at the maximum operational limit).
- `resources/<run_name>/additional_heating_indicator.nc`
    NetCDF file containing a binary indicator for additional after heating requirements.
"""

import logging
import xarray as xr

logger = logging.getLogger(__name__)

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from _helpers import set_scenario_config
from tes_additional_heating_approximator import PTESAdditionalHeatingApproximator
from build_TES_temperature_profiles import BuildTESTemperature

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake(
            "build_ptes_temperature_profiles",
            clusters=6,
            planning_horizons="2030",
        )

    set_scenario_config(snakemake)

    # Load forward temperature profile
    logger.info("Loading forward temperature profile")
    forward_temp = xr.open_dataarray(snakemake.input.central_heating_forward_temperature_profiles)

    # Load parameter for max. PTES temperature.
    max_PTES_temperature = snakemake.params.max_PTES_temperature
    logger.info(f"Using maximum PTES direct usage temperature: {max_PTES_temperature}°C")

    # Initialize TES top temperature approximator.
    tes_temp_builder = BuildTESTemperature(
        forward_temperature_celsius=forward_temp,
        max_PTES_temperature=max_PTES_temperature,
    )

    # Approximate PTES top temperature profile
    ptes_top_temperature = tes_temp_builder.clipped_top_temperature

    # Initialize additional heating approximator.
    heating_approximator = PTESAdditionalHeatingApproximator(
        forward_temperature_celsius=forward_temp,
        max_PTES_temperature=max_PTES_temperature,
    )

    additional_heating = heating_approximator.determine_ptes_usage()

    # Save output
    logger.info(f"Saving PTES top temperature profile to {snakemake.output.ptes_top_temperature}")
    ptes_top_temperature.to_netcdf(snakemake.output.ptes_top_temperature)

    # Save output
    logger.info(f"Saving additional heating indicator to {snakemake.output.additional_heating_indicator}")
    additional_heating.to_netcdf(snakemake.output.additional_heating_indicator)