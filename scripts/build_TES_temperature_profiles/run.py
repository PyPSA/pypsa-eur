# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Determine PTES top temperature and additional heating requirements for integration.

This script evaluates the district heating forward temperature profile to:
  1. Approximate the top temperature of the Pit Thermal Energy Storage (PTES) system,
     ensuring that the temperature does not exceed the operational limit.
  2. Determine whether additional (after) heating is needed. A binary indicator is generated:
       - 0: The forward temperature is less than or equal to the PTES maximum; direct usage is possible.
       - 1: The forward temperature exceeds the PTES maximum; additional heating (e.g., via a heat pump) is required.

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
import sys
import os

# Set up logging.
logger = logging.getLogger(__name__)

# Ensure the project root directory is on the Python path.
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

    # Set scenario-specific configuration.
    set_scenario_config(snakemake)

    # Load the forward temperature profile for district heating.
    logger.info("Loading forward temperature profile")
    forward_temp = xr.open_dataarray(snakemake.input.central_heating_forward_temperature_profiles)

    # Load the configuration parameter for PTES maximum direct usage temperature.
    max_PTES_temperature = snakemake.params.max_PTES_temperature
    logger.info(f"Using maximum PTES direct usage temperature: {max_PTES_temperature}°C")

    # Initialize the TES top temperature approximator.
    tes_temp_builder = BuildTESTemperature(
        forward_temperature_celsius=forward_temp,
        max_PTES_temperature=max_PTES_temperature,
    )
    # Approximate the PTES top temperature profile (clipped at max temperature).
    ptes_top_temperature = tes_temp_builder.clipped_top_temperature

    # Initialize the additional heating approximator.
    heating_approximator = PTESAdditionalHeatingApproximator(
        forward_temperature_celsius=forward_temp,
        max_PTES_temperature=max_PTES_temperature,
    )
    # Determine the additional heating indicator.
    # 0 indicates direct use of stored heat; 1 indicates additional after heating is required.
    additional_heating = heating_approximator.determine_ptes_usage()

    # Save the PTES top temperature profile to its own NetCDF file.
    logger.info(f"Saving PTES top temperature profile to {snakemake.output.ptes_top_temperature}")
    ptes_top_temperature.to_netcdf(snakemake.output.ptes_top_temperature)

    # Save the additional heating indicator to its own NetCDF file.
    logger.info(f"Saving additional heating indicator to {snakemake.output.additional_heating_indicator}")
    additional_heating.to_netcdf(snakemake.output.additional_heating_indicator)





# was wir als output brauchen ist ein Temperautr profil für PTES welches wir verwenden, welches aber bei der konstaten max temp liegt
# was wir dann brauchen ist eine Möglichkeit um zu integrieren dass wir mit PTES max Temp arbeiten können
# temperatur des Speichers ist T_forward und T_return, wobei wenn T_forward > PTES_max dann ist PTES_max die T_top temp und wir brauchen eine Nacherhitzung
# es brauch eine Funktion für den COP, bzw eine Möglichkeit dass eine Quelle zu gewissen Zeiten nicht für die WP verfügbar sind