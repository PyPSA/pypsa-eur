# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Determine additional heating requirements for PTES integration.

This script evaluates whether additional (after) heating is needed based on a comparison
between the district heating forward temperature profile and the maximum direct usage
temperature of Pit Thermal Energy Storage (PTES). A binary indicator is generated:
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
    Forward temperature profile (in Celsius).

Outputs
-------
- `resources/<run_name>/additional_heating_indicator.nc`
    NetCDF file containing a binary indicator for additional after heating requirements.
"""

import logging
import xarray as xr
import sys
import os

# Set up logging
logger = logging.getLogger(__name__)

# Ensure the project root directory is on the Python path.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from _helpers import set_scenario_config
from tes_additional_heating_approximator import PTESAdditionalHeatingApproximator

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake(
            "build_ptes_additional_heating_indicator",
            clusters=6,
            planning_horizons="2030",
        )

    # Set scenario-specific configuration.
    set_scenario_config(snakemake)

    # Load the forward temperature profile for district heating.
    logger.info("Loading forward temperature profile")
    forward_temp = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )

    # Load the configuration parameter for PTES maximum direct usage temperature.
    max_PTES_temperature = snakemake.params.max_PTES_temperature
    logger.info(f"Using maximum PTES direct usage temperature: {max_PTES_temperature}°C")

    # Initialize the additional heating approximator.
    approximator = PTESAdditionalHeatingApproximator(
        forward_temperature_celsius=forward_temp,
        max_PTES_temperature=max_PTES_temperature,
    )

    # Determine the additional heating indicator.
    # 0 indicates the stored heat meets the demand directly;
    # 1 indicates that after heating is required.
    additional_heating = approximator.determine_ptes_usage()

    logger.info(f"Saving additional heating indicator to {snakemake.output.additional_heating_indicator}")
    # Save the binary indicator to a NetCDF file.
    additional_heating.to_netcdf(snakemake.output.additional_heating_indicator)




# was wir als output brauchen ist ein Temperautr profil für PTES welches wir verwenden, welches aber bei der konstaten max temp liegt
# was wir dann brauchen ist eine Möglichkeit um zu integrieren dass wir mit PTES max Temp arbeiten können
# temperatur des Speichers ist T_forward und T_return, wobei wenn T_forward > PTES_max dann ist PTES_max die T_top temp und wir brauchen eine Nacherhitzung
# es brauch eine Funktion für den COP, bzw eine Möglichkeit dass eine Quelle zu gewissen Zeiten nicht für die WP verfügbar sind