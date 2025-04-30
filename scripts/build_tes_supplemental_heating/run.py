# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Approximate thermal energy storage (TES) top temperature and identify need for supplemental heating.

This script evaluates the district heating forward temperature profile to:
  - Approximate the top temperature of the TES system,
     ensuring that the temperature does not exceed the operational limit.
  - Determine whether supplemental heating is needed. A binary indicator is generated:
       - 1: The forward temperature is less than or equal to the TES maximum; direct usage is possible.
       - 0: The forward temperature exceeds the TES maximum; supplemental heating (e.g., via a heat pump) is required.

Relevant Settings
-----------------
.. code:: yaml
    sector
        district_heating:
            ptes:
                supplemental_heating:
                    enable: true
                max_top_temperature: 90  # Maximum PTES temperature (°C) usable without supplemental heating

Inputs
------
- `resources/<run_name>/forward_temperature.nc`
    Forward temperature profile for the district heating network.

Outputs
-------
- `resources/<run_name>/tes_top_temperature_profile.nc`
    Clipped TES top temperature profile (in °C).
- `resources/<run_name>/additional_heating_indicator.nc`
    Binary indicator for additional heating (1 = direct TES use, 0 = supplemental heating required).
"""

import logging
import xarray as xr

logger = logging.getLogger(__name__)

from _helpers import set_scenario_config
from tes_supplemental_heating_approximator import TESSupplementalHeatingApproximator
from build_tes_top_temperature_profile import BuildTesTopTemperature

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
    tes_temp_builder = BuildTesTopTemperature(
        forward_temperature_celsius=forward_temp,
        max_PTES_temperature=max_PTES_temperature,
    )

    # Approximate PTES top temperature profile
    ptes_top_temperature = tes_temp_builder.clipped_top_temperature

    # Initialize supplemental heating approximator.
    supplemental_heating_approximator = TESSupplementalHeatingApproximator(
        forward_temperature_celsius=forward_temp,
        max_PTES_temperature=max_PTES_temperature,
    )

    supplemental_heating_profile = supplemental_heating_approximator.determine_ptes_usage()

    # Save output
    logger.info(f"Saving TES top temperature profile to {snakemake.output.tes_top_temperature_profile}")
    ptes_top_temperature.to_netcdf(snakemake.output.tes_top_temperature_profile)

    # Save output
    logger.info(f"Saving supplemental heating profile to {snakemake.output.tes_supplemental_heating_profile}")
    supplemental_heating_profile.to_netcdf(snakemake.output.tes_supplemental_heating_profile)
